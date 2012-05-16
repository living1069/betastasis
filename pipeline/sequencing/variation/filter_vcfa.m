
function [] = filter_vcfa(vcfa_file, out_file, varargin)

%min_reads = 5;
discard_silent = false;
%strand_bias_threshold = 0.001;
abs_negative_groups = [];
sample_map = [];

groups = {};

% Use a scoring system based on a ranksum test on the phred-scaled
% genotype likelihoods.
ranksum_groups = [];
ranksum_min_significance = 0.50;     % P-value lower or variant is discarded

cosmic_score = 2;
kgenomes_score = -1;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'sample.*map')
		sample_map = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'cosmic.*score')
		cosmic_score = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'kgenome.*score')
		kgenomes_score = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'ranksum.*group')
		ranksum_groups = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'ranksum.*signif')
		ranksum_min_significance = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'group')
		groups = varargin{k+1};
		% Odd cells are group names, even cells are index/logical.
		for g = 1:2:length(groups)
			if ~ischar(groups{g}), error 'Invalid group format.'; end
			if islogical(groups{g+1})
				groups{g+1} = find(groups{g+1});
			elseif ~isnumeric(groups{g+1})
				error 'Invalid group format.';
			end
		end
		continue;
	end
	
	%if regexpi(varargin{k}, 'min.*reads')
%		min_reads = varargin{k+1}; continue;
%	end
	
	if regexpi(varargin{k}, 'discard.*silent')
		discard_silent = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'abs.*neg')
		abs_negative_groups = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

[data, headers] = readtable(vcfa_file);
if any(rx(headers, 'HOM/HET'))
	error 'File has already been filtered!';
end

first_sample_col = find(rx(headers, 'COSMIC'))+1;
S = length(headers) - first_sample_col + 1;
V = length(data{1});

samples = headers(first_sample_col:end);
%samples = regexprep(samples, '.*/', '');
%samples = regexprep(samples, '\.bam$', '');

% Relabel the samples.
if ~isempty(sample_map)
	if strcmp(class(sample_map), 'containers.Map')
		valid = sample_map.isKey(samples);
		samples(valid) = sample_map.values(samples(valid));
	elseif iscellstr(sample_map) && length(sample_map) == length(samples)
		samples = sample_map;
	else
		error 'Invalid sample map format.';
	end
end

% Permute the samples if the user has specified groups.
if ~isempty(groups)
	sample_order = [];
	
	for g = 1:2:length(groups)
		gsamples = groups{g+1};
		gsamples = gsamples(:)';
		groups{g+1} = (1:length(gsamples))+length(sample_order);
		sample_order = [sample_order gsamples];
	end
	if length(sample_order) ~= length(samples)
		error 'Specified groups do not cover all samples.';
	end
	samples = samples(sample_order);
	headers(first_sample_col:end) = samples;
	data(first_sample_col:end) = data(sample_order+first_sample_col-1);
end





% Construct numeric matrices out of the variant-sample spreadsheet.
consequences = data{rx(headers, 'SYNONYMOUS')};
total_reads = nan(V, S);
genotype = nan(V, S);
kgenomes = ~strcmpi(data{rx(headers, '1000GENOMES')}, '-');
cosmic = ~strcmpi(data{rx(headers, 'COSMIC')}, '-');

for s = 1:S
	info = data{first_sample_col + s - 1};
	tokens = regexp(info, '(.+?):.+?:(\d+)', 'tokens');
	for v = 1:length(tokens)
		token = tokens{v}{1};
		if strcmp(token{1}, './.')
			genotype(v, s) = NaN;
		else
			genotype(v, s) = sum(token{1} == '1');
		end
		total_reads(v, s) = str2double(token{2});
	end
end





keep = true(V, 1);

% Filter out silent variants.
if discard_silent
	keep = keep & ~rx(consequences, '^(synonymous|unknown|-)');
end

if ~isempty(abs_negative_groups)
	abs_negative_samples = unique([groups{abs_negative_groups*2}])
	keep = keep & ~any(genotype(:, abs_negative_samples) > 0, 2);
end

if ~isempty(ranksum_groups)
	% Score the variants by t-testing for difference between test and
	% reference samples.
	test_samples = unique([groups{ranksum_groups{1}*2}]);
	ref_samples = unique([groups{ranksum_groups{2}*2}]);
	
	mutated = genotype;
	mutated(mutated == 2) = 1;
	
	[~, p] = ttest2(mutated(:, test_samples)', mutated(:, ref_samples)', ...
		0.05, 'right');
	
	score = (1 - p)';
	keep = keep & (p <= ranksum_min_significance)';

else
	% Just score the variants based on frequency.
	score = nansum(genotype / 2, 2) / size(genotype, 2);
end



% Adjust the scoring based on whether the variant is reported by the COSMIC
% and 1000 Genomes projects.
score(cosmic) = score(cosmic) + cosmic_score;
score(kgenomes) = score(kgenomes) + kgenomes_score;



[~, order] = sort(score, 'descend');
order = order(keep(order));








% Figure out where we need to insert the homozygous/heterozygous counts.
fid = fopen(out_file, 'W');
fprintf(fid, '%s\t', headers{1:(first_sample_col-1)});

% Print headers for heterozygous/homozygous totals
if isempty(groups)
	fprintf(fid, 'HOM/HET');
else
	for g = 1:2:length(groups)
		fprintf(fid, '%s HOM/HET\t', upper(groups{g}));
	end
end

fprintf(fid, '%s\t', headers{first_sample_col:end-1});
fprintf(fid, '%s\n', headers{end});

genotype_strs = { '0/0', '0/1', '1/1' };

for k = order'
	for c = 1:(first_sample_col-1)
		fprintf(fid, '%s\t', data{c}{k});
	end
	
	if isempty(groups)
		fprintf(fid, '%d / %d', sum(genotype(k, :) == 2), ...
			sum(genotype(k, :) > 0));
	else
		for g = 1:2:length(groups)
			fprintf(fid, '%d / %d', sum(genotype(k, groups{g+1}) == 2), ...
				sum(genotype(k, groups{g+1}) > 0));
			if g < length(groups) - 1, fprintf(fid, '\t'); end
		end
	end
		
	for s = 1:size(genotype, 2)
		if ~isnan(genotype(k, s))
			fprintf(fid, '\t%s (%d)', genotype_strs{genotype(k, s)+1}, ...
				total_reads(k, s));
		else
			fprintf(fid, '\t');
		end
	end
	fprintf(fid, '\n');
end
fclose(fid);


