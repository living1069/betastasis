function [] = copa_analysis(expr, qtile)

global organism;
genes = organism.Genes;

norm_expr = log2(expr.Mean);
norm_expr = norm_expr - repmat(median(norm_expr, 2), 1, size(norm_expr, 2));
norm_expr = norm_expr ./ ...
	repmat(median(abs(norm_expr), 2), 1, size(norm_expr, 2));

qtile_expr = quantile(norm_expr, qtile, 2);

valid = find(~isnan(qtile_expr));

% Identify the samples with the strongest XYth percentile expression.
if qtile > 0.5
	[~, order] = sort(qtile_expr(valid), 'descend');
elseif qtile < 0.5
	[~, order] = sort(qtile_expr(valid), 'ascend');
end

order = valid(order);

fprintf(1, 'Genes with strongest %.2f quantile expression:\n', qtile);
for g = order(1:20)'
	fprintf(1, '%s (%f)\n', genes.Name{g}, qtile_expr(g));
end



