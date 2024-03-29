function fmatrix = fmatrix_diff_expression(test, ref, varargin)

global organism;
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
fmatrix.features = strcat('N:EXPR:', expr.rows.gene);
fmatrix.data = log2(expr.mean);

fprintf('%d differential gene expression features.\n', F);

