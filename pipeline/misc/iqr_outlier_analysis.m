function [] = iqr_outlier_analysis(expr)

global organism;
genes = organism.Genes;

log_expr = log2(expr.Mean);
S = size(log_expr, 2);
G = size(log_expr, 1);

quartiles = quantile(log_expr, [0.25 0.75], 2);
q25 = quartiles(:, 1);
q75 = quartiles(:, 2);
iqr = q75 - q25;

high_outliers = nan(G, 1);
low_outliers = nan(G, 1);

for g = 1:length(genes.Name)
	high_outliers(g) = sum(log_expr(g, :) > q75(g) + 1.5 * iqr(g));
	low_outliers(g) = sum(log_expr(g, :) < q25(g) - 1.5 * iqr(g));
end

valid = find(any(~isnan(log_expr), 2));

score = high_outliers + low_outliers;
whos score

% Identify the genes with the largest number of samples with outlier expression.
[~, order] = sort(score(valid), 'descend');
order = valid(order);

fprintf(1, 'Genes with largest number of outliers (IQR):\n');
for g = order(1:20)'
	fprintf(1, '%s (%.2f%%)\n', genes.Name{g}, score(g) / S * 100);
end



