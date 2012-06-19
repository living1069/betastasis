function fmatrix = fmatrix_gene_expression(expr, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = expr.meta.sample_id;
if length(samples) ~= length(unique(samples))
	error 'Dataset contains technical replicates.';
end

valid = find(~any(isnan(expr.mean), 2));
expr = filter_rows(expr, valid);

S = length(expr.meta.sample_id);
F = length(valid);

fmatrix.samples = expr.meta.sample_id;
fmatrix.features = strcat('N:EXPR:', expr.rows.gene_symbol);
fmatrix.data = log2(expr.mean);

fprintf('%d gene expression features.\n', F);

