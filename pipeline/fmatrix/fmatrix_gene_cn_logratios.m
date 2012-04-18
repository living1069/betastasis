function fmatrix = fmatrix_gene_cn_logratios(test, ref, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = test.meta.sample_id;
if length(samples) ~= length(unique(samples))
	error 'Dataset contains technical replicates.';
	%fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	%segments = merge_replicates(segments);
end

gene_lr = gene_cgh_logratios(test, ref, 'MinRegion', 100e3);

valid = find(~any(isnan(gene_lr), 2));

S = length(test.meta.sample_id);
F = length(valid);

fmatrix.Samples = samples;
fmatrix.Features = cell(F, 1);
fmatrix.Data = nan(F, S);

for f = 1:length(valid)
	g = valid(f);
	
	if isnan(genes.Chromosome(g))
		pos_str = '?:?:?:?';
	else
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{genes.Chromosome(g)}, ...
			genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
	end

	fmatrix.Features{f} = sprintf('N:CNLR:%s:%s', genes.Name{g}, pos_str);
	fmatrix.Data(f, :) = gene_lr(g, :);
end

fprintf(1, '%d copy number logratio features.\n', F);


