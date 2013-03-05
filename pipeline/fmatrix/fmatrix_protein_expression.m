function fmatrix = fmatrix_protein_expression(expr, varargin)

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
fmatrix.features = strcat('N:EXPR:', expr.rows.protein);
fmatrix.data = expr.mean;

fprintf('%d protein expression features.\n', F);

