function top = top_fold_change(expr, ref_expr, features, ...
	ratio_threshold, abs_threshold)

% The log-2 ratio threshold defaults to 1. A feature is considered
% differentially expressed if the log ratio between the reference and
% measured expression is outside the range [-threshold, threshold].
if nargin < 4
	ratio_threshold = 1;
end

if nargin < 5
	abs_threshold = 0;
end

if isfield(expr, 'Mean'), expr = expr.Mean; end
if isfield(ref_expr, 'Mean'), ref_expr = ref_expr.Mean; end

expr = median(expr, 2);
ref_expr = median(ref_expr, 2);

log_ratios = log2(expr ./ ref_expr);
sum_expr = expr + ref_expr;

positive_diff_expr = find((log_ratios > ratio_threshold) & ...
	(sum_expr > abs_threshold) & (log_ratios ~= Inf));
negative_diff_expr = find((log_ratios < -ratio_threshold) & ...
	(sum_expr > abs_threshold) & (log_ratios ~= -Inf));

[~, order] = sort(log_ratios(positive_diff_expr), 'descend');
positive_diff_expr = positive_diff_expr(order);

[~, order] = sort(log_ratios(negative_diff_expr), 'ascend');
negative_diff_expr = negative_diff_expr(order);


fprintf(1, 'Positively differentially expressed features:\n');
for k = 1:length(positive_diff_expr)
	fprintf(1, '- %s, log ratio %.2f (%.2f vs %.2f)\n', ...
		features.Name{positive_diff_expr(k)}, ...
		log_ratios(positive_diff_expr(k)), ...
		expr(positive_diff_expr(k)), ...
		ref_expr(positive_diff_expr(k)));
end

fprintf(1, 'Negatively differentially expressed features:\n');
for k = 1:length(negative_diff_expr)
	fprintf(1, '- %s, log ratio %.2f (%.2f vs %.2f)\n', ...
		features.Name{negative_diff_expr(k)}, ...
		log_ratios(negative_diff_expr(k)), ...
		expr(negative_diff_expr(k)), ...
		ref_expr(negative_diff_expr(k)));
end

top = struct;
top.Positive = positive_diff_expr;
top.Negative = negative_diff_expr;
top.LogRatio = log_ratios;

