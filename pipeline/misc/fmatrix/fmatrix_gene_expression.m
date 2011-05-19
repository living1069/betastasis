function fmatrix = fmatrix_gene_expression(expr, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = expr.Meta.Sample.ID;
if length(samples) ~= length(unique(samples))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	expr = merge_replicates(expr);
end

valid = find(~any(isnan(expr.Mean), 2));

S = length(expr.Meta.Sample.ID);
F = length(valid);

fmatrix.Samples = expr.Meta.Sample.ID;
fmatrix.Features = cell(F, 1);
fmatrix.Data = log2(expr.Mean(valid, :));

for f = 1:F
	g = valid(f);
	
	if isnan(genes.Chromosome(g))
		pos_str = '?:?:?:?';
	else
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{genes.Chromosome(g)}, ...
			genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
	end
	
	fmatrix.Features{f} = sprintf('N:EXPR:%s:%s', genes.Name{g}, pos_str);
end

fprintf(1, '%d gene expression features.\n', F);

