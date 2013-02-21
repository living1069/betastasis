
function norm = normalize_ratio_medians(expr, detect_threshold)

norm = expr;

valid = median(expr.mean, 2) > detect_threshold;

fprintf('Using %d features for median-of-ratios normalization.\n', sum(valid));

ratios = nanmedian(expr.mean(valid, :) ./ ...
	repmat(expr.mean(valid, 1), 1, size(expr.mean, 2)), 1);
ratios = ratios / nanmedian(ratios);

for s = 1:size(norm.mean, 2)
	norm.mean(:, s) = norm.mean(:, s) / ratios(s);
end

