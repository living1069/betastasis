
function [] = filter_vcfa(vcfa_file, out_file, varargin)

global organism;

%strand_bias_threshold = 0.001;
min_genotype_qual = 30;
abs_negative_groups = [];
sample_map = [];

groups = {};

% Use a scoring system based on a ranksum test on the phred-scaled
% genotype likelihoods.
ranksum_groups = [];
ranksum_min_significance = 0.50;     % P-value lower or variant is discarded

cosmic_score = 2;
kgenomes_score = -1;
silent_score = -Inf;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'sample.*map')
		sample_map = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'geno.*qual')
		min_genotype_qual = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'silent.*score')
		silent_score = varargin{k+1}; continue;
	end
	if rx(varargin{k}, 'cosmic.*score')
		cosmic_score = varargin{k+1}; continue;
	end
	if rx(varargin{k}, 'kgenome.*score')
		kgenomes_score = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'abs.*neg')
		abs_negative_groups = varargin{k+1}; continue;
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
	
	error('Unrecognized option "%s".', varargin{k});
end

variants = read_vcf(vcfa_file)

S = length(variants.meta.sample_id);
V = size(variants.genotype, 1);

% Relabel the samples.
if ~isempty(sample_map)
	if strcmp(class(sample_map), 'containers.Map')
		valid = sample_map.isKey(variants.meta.sample_id);
		variants.meta.sample_id(valid) = ...
			sample_map.values(variants.meta.sample_id(valid));
	elseif iscellstr(sample_map) && length(sample_map) == S
		variants.meta.sample_id = sample_map;
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
	if length(sample_order) ~= S
		error 'User specified groups do not cover all samples.';
	end
	variants = filter(variants, sample_order);
end





% Set genotypes to unknown if their quality is not sufficient.
variants.genotype(~(variants.genotype_quality >= min_genotype_qual)) = NaN;
silent = rx(variants.rows.synonymous, '^(synonymous|unknown|-)');



keep = true(V, 1);

% Filter out silent variants if bonus for coding variants is infinite.
if silent_score == -Inf
	keep = keep & ~silent;
end

if ~isempty(abs_negative_groups)
	abs_negative_samples = unique([groups{abs_negative_groups*2}])
	keep = keep & ~any(variants.genotype(:, abs_negative_samples) > 0, 2);
end

if ~isempty(ranksum_groups)
	% Score the variants by t-testing for difference between test and
	% reference samples.
	test_samples = unique([groups{ranksum_groups{1}*2}]);
	ref_samples = unique([groups{ranksum_groups{2}*2}]);
	
	mutated = variants.genotype;
	mutated(mutated == 2) = 1;
	
	[~, p] = ttest2(mutated(:, test_samples)', mutated(:, ref_samples)', ...
		0.05, 'right');
	
	score = (1 - p)';
	keep = keep & (p <= ranksum_min_significance)';

else
	% Just score the variants based on frequency.
	score = nansum(variants.genotype / 2, 2) / size(variants.genotype, 2);
end



% Adjust the scoring based on whether the variant is reported by the COSMIC
% and 1000 Genomes projects.
if isfinite(silent_score)
	score = score + silent_score * silent;
end
score = score + cosmic_score * variants.rows.cosmic;
score = score + kgenomes_score * ~strcmp(variants.rows.kgenomes, '-');



[~, order] = sort(score, 'descend');
order = order(keep(order));








% Figure out where we need to insert the homozygous/heterozygous counts.
fid = fopen(out_file, 'W');
fprintf(fid, ['CHROMOSOME\tPOSITION\tREFERENCE\tOBSERVED\tFUNCTION\t' ...
	'NEARBY_GENES\tSYNONYMOUS\tPROTEIN_EFFECT\tDBSNP\t1000GENOMES\tCOSMIC']);

% Print headers for heterozygous/homozygous totals
if isempty(groups)
	fprintf(fid, '\tHOM/HET/TOT');
else
	for g = 1:2:length(groups)
		fprintf(fid, '\t%s HOM/HET/TOT', upper(groups{g}));
	end
end

fprintf(fid, '\t%s', variants.meta.sample_id{:});
fprintf(fid, '\n');

chr_names = organism.Chromosomes.Name;
r = variants.rows;

genotype_str = { '0/0', '0/1', '1/1' };

kgenome_str = {''; 'YES'};
kgenome_str = kgenome_str((r.kgenomes >= 0) + 1);

cosmic_str = {''; 'YES'};
cosmic_str = cosmic_str(r.cosmic + 1);

for k = order'
	fprintf(fid, 'chr%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t', ...
		chr_names{r.chromosome(k)}, r.position(k), r.ref_allele{k}, ...
		r.alt_allele{k}, r.function{k}, r.nearby_genes{k}, r.synonymous{k}, ...
		r.protein_effect{k}, r.dbsnp{k}, kgenome_str{k}, cosmic_str{k});

	if isempty(groups)
		fprintf(fid, '%d / %d / %d', sum(variants.genotype(k, :) == 2), ...
			sum(variants.genotype(k, :) > 0), ...
			sum(~isnan(variants.genotype(k, :))));
	else
		for g = 1:2:length(groups)
			fprintf(fid, '%d / %d / %d', ...
				sum(variants.genotype(k, groups{g+1}) == 2), ...
				sum(variants.genotype(k, groups{g+1}) > 0), ...
				sum(~isnan(variants.genotype(k, groups{g+1}))));
			if g < length(groups) - 1, fprintf(fid, '\t'); end
		end
	end
	
	for s = 1:S
		if ~isnan(variants.genotype(k, s))
			fprintf(fid, '\t%s (%d)', ...
				genotype_str{variants.genotype(k,s)+1}, ...
				variants.genotype_quality(k, s));
		else
			fprintf(fid, '\t');
		end
	end
	fprintf(fid, '\n');
end
fclose(fid);


