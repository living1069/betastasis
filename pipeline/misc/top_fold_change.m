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
sum_expr = (expr + ref_expr) / 2;

diff_expr = find((abs(log_ratios) > ratio_threshold) & ...
	(sum_expr > abs_threshold) & (abs(log_ratios) ~= Inf));

[~, order] = sort(abs(log_ratios(diff_expr)), 'descend');
diff_expr = diff_expr(order);

fprintf(1, 'Feature\tLogratio\tAbs logratio\tTest expr\tRef expr\n');
for k = diff_expr'
	fprintf(1, '%s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
		features.Name{k}, log_ratios(k), abs(log_ratios(k)), ...
		expr(k), ref_expr(k));
end

