function [] = top_expr(expr, features)

if isfield(expr, 'Mean'), expr = expr.Mean; end
expr = median(expr, 2);

notnan = find(~isnan(expr));
[~, order] = sort(expr, 'descend');

fprintf(1, 'Most highly expressed genes:\n');
for k = 1:20
	idx = notnan(order(k));
	fprintf(1, '- %s, expr %.2f, log2 expr %.2f\n', ...
		features.Name{idx}, expr(idx), log2(expr(idx)));
end

