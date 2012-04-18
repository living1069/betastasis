
function [] = variant_overlap(vcfa_file, varargin)

min_reads = 5;
discard_silent = false;
strand_bias_threshold = 0.001;

groups = {};

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'group')
		groups = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'min.*reads')
		min_reads = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'discard.*silent')
		discard_silent = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'abs.*neg')
		abs_negative_groups = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

[data, headers] = readtable(vcfa_file);
if any(rx(headers, 'HET/HOM'))
	error 'File has already been filtered!';
end

first_sample_col = find(rx(headers, 'VARIANT_CALL_DETAILS'))+1;
S = length(headers) - first_sample_col + 1;
V = length(data{1});

samples = headers(first_sample_col:end);
samples = regexprep(samples, '.*/', '');
samples = regexprep(samples, '\.bam$', '');

% Relabel the samples.
if ~isempty(sample_map)
	valid = sample_map.isKey(samples);
	samples(valid) = sample_map.values(samples(valid));
end

% Permute the samples if the user has specified groups.
if ~isempty(groups)
	group_samples = {};
	sample_order = [];
	for g = 1:length(groups)
		gsamples = find(rx(samples, groups{g}));
		group_samples{g} = (1:length(gsamples))+length(sample_order);
		sample_order = [sample_order gsamples'];
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
strand_bias = nan(V, 1);
kgenomes = ~strcmpi(data{rx(headers, '1000GENOMES')}, '-');
cosmic = ~strcmpi(data{rx(headers, 'COSMIC')}, '-');

% Parse strand bias.
details = data{rx(headers, 'VARIANT_CALL_DETAILS')};
tokens = regexp(details, 'PV4=(.+?),', 'tokens');
for v = 1:length(tokens)
	if isempty(tokens{v}), strand_bias(v) = 1; continue; end
	strand_bias(v) = str2double(tokens{v}{1}{1});
end


for s = 1:S
	info = data{first_sample_col + s - 1};
	tokens = regexp(info, '(.+?):.+?:(\d+):.*', 'tokens');
	for v = 1:length(tokens)
		token = tokens{v}{1};
		genotype(v, s) = sum(token{1} == '1');
		total_reads(v, s) = str2double(token{2});
	end
end






keep = true(V, 1);

% Filter out genotypes that are based on less than the minimum amount of reads.
if ~isempty(min_reads)
	genotype(total_reads < min_reads) = NaN;
	keep = keep & any(genotype > 0, 2);
end

% Filter out silent variants.
if discard_silent
	keep = keep & ~rx(consequences, '^(synonymous|unknown|-)');
end

% Filter out variants with strand bias.
if strand_bias_threshold > 0
	keep = keep & (strand_bias >= strand_bias_threshold);
end

if ~isempty(abs_negative_groups)
	abs_negative_samples = unique([group_samples{abs_negative_groups}])
	keep = keep & ~any(genotype(:, abs_negative_samples) > 0, 2);
end

if ~isempty(ranksum_groups)
	% Score the variants by t-testing for difference between test and
	% reference samples.
	test_samples = unique([group_samples{ranksum_groups{1}}]);
	ref_samples = unique([group_samples{ranksum_groups{2}}]);
	
	[~, p] = ttest2(genotype(:, test_samples)', genotype(:, ref_samples)', ...
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
counts_after_col = find(rx(headers, 'VARIANT_CALL_DETAILS'));

fid = fopen(out_file, 'W');
fprintf(fid, '%s\t', headers{1:counts_after_col});
fprintf(fid, 'HET/HOM\tHOM\t');
fprintf(fid, '%s\t', headers{counts_after_col+1:end-1});
fprintf(fid, '%s\n', headers{end});

genotype_strs = { '0/0', '0/1', '1/1' };

for k = order'
	for c = 1:counts_after_col
		fprintf(fid, '%s\t', data{c}{k});
	end
	fprintf(fid, '%d\t%d', sum(genotype(k, :) > 0), sum(genotype(k, :) == 2));
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


