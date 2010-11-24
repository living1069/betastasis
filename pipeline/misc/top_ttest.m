function p = top_ttest(expr, ref, features, p_threshold, abs_threshold)

if nargin < 4
	p_threshold = 0.05;
end
if nargin < 5
	abs_threshold = 0;
end

if isfield(expr, 'Mean'), expr = expr.Mean; end
if isfield(ref, 'Mean'), ref = ref.Mean; end

whos expr ref;	

[~, p] = ttest2(expr', ref');

expr_mean = nanmean(expr, 2);
ref_mean = nanmean(ref, 2);

sig = find((p <= p_threshold)' & ...
	max([expr_mean ref_mean], [], 2) > abs_threshold);

[~, order] = sort(p(sig), 'ascend');
sig = sig(order);

fprintf(1, 'Differentially expressed features:\n');
for k = 1:length(sig)
	fprintf(1, '- %s, p-value %g (%.2f vs %.2f)\n', ...
		features.Name{sig(k)}, p(sig(k)), ...
		expr_mean(sig(k)), ref_mean(sig(k)));
end

