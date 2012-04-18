function fmatrix = fmatrix_gene_expression(expr, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = expr.meta.sample_id;
if length(samples) ~= length(unique(samples))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	expr = merge_replicates(expr);
end

valid = find(~any(isnan(expr.mean), 2));

S = length(expr.meta.sample_id);
F = length(valid);

fmatrix.samples = expr.meta.sample_id;
fmatrix.features = cell(F, 1);
fmatrix.data = log2(expr.mean(valid, :));

for f = 1:F
	g = valid(f);
	
	if isnan(genes.Chromosome(g))
		pos_str = '?:?:?:?';
	else
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{genes.Chromosome(g)}, ...
			genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
	end
	
	fmatrix.features{f} = sprintf('N:EXPR:%s:%s', genes.Name{g}, pos_str);
end

fprintf(1, '%d gene expression features.\n', F);

