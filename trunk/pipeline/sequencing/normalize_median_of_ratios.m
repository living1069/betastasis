
function [norm, ratios] = normalize_median_of_ratios(expr, detect_threshold)

norm = nan(size(expr));

valid = median(expr, 2) > detect_threshold;

fprintf('Using %d features for median-of-ratios normalization.\n', sum(valid));

ratios = nanmedian(expr(valid,:) ./ repmat(expr(valid,1), 1, size(expr,2)), 1);
ratios = ratios / nanmedian(ratios);
for s = 1:size(norm, 2), norm(:, s) = expr(:, s) / ratios(s); end

