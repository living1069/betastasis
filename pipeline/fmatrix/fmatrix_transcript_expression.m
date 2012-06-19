function fmatrix = fmatrix_gene_expression(expr, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = expr.meta.sample_id;
if length(samples) ~= length(unique(samples))
	error 'Dataset contains technical replicates.';
end

valid = find(~any(isnan(expr.mean), 2));

S = length(expr.meta.sample_id);
F = length(valid);

fmatrix.samples = expr.meta.sample_id;
fmatrix.features = strcat('N:EXPR:', expr.rows.gene_symbol, ...
	repmat({' ('}, F, 1), expr.rows.transcript, repmat({'):'}, F, 1));
fmatrix.data = log2(expr.mean(valid, :));

fprintf('%d transcript expression features.\n', F);

