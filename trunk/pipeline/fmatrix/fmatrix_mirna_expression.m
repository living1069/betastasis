function fmatrix = fmatrix_mirna_expression(expr, varargin)

samples = expr.meta.sample_id;
if length(samples) ~= length(unique(samples))
	error('Dataset contains technical replicates.');
end

valid = find(any(~isnan(expr.mean), 2));

S = length(expr.meta.sample_id);
F = length(valid);

fmatrix.samples = expr.meta.sample_id;
fmatrix.features = cell(F, 1);
fmatrix.data = log2(expr.mean(valid, :));

for f = 1:F
	g = valid(f);
	fmatrix.features{f} = sprintf('N:EXPR:%s', expr.rows.mirna{g});
end

fprintf('%d microRNA expression features.\n', F);

