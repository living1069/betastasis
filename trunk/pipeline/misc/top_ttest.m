function p = top_ttest(expr, ref, features, p_threshold, abs_threshold)

if nargin < 4
	p_threshold = 0.05;
end
if nargin < 5
	abs_threshold = 0;
end

if isfield(expr, 'Mean'), expr = expr.Mean; end
if isfield(ref, 'Mean'), ref = ref.Mean; end

if any(any(expr <= 0)) || any(any(ref <= 0))
	fprintf(1, 'Non-positive expression values found. Adding background...\n');
	expr = expr - min(min(expr)) + 1;
	ref = ref - min(min(ref)) + 1;
end
	
log_expr = log2(expr);
log_ref = log2(ref);

[~, p] = ttest2(log_expr', log_ref');

expr_mean = nanmean(log_expr, 2);
ref_mean = nanmean(log_ref, 2);

sig = find((p <= p_threshold)' & ...
	max([expr_mean ref_mean], [], 2) > abs_threshold);

[~, order] = sort(p(sig), 'ascend');
sig = sig(order);

fprintf(1, 'Feature\tp-value\tTest expr\tRef expr\n');
for k = sig'
	fprintf(1, '%s\t%g\t%.2f\t%.2f\n', ...
		features.Name{k}, p(k), expr_mean(k), ref_mean(k));
end

