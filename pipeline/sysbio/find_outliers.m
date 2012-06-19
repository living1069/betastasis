
% Author: Matti Annala <matti.annala@tut.fi>

function aberrated = find_outliers(test, ref, ...
	logratio_thresholds, abs_thresholds)

ref_percentile = 1.0;

if nargin < 3, logratio_thresholds = [-3 -2 -1 1 2 3]; end
if nargin < 4, abs_thresholds = [0 0 0 0 0 0]; end
	
if any(logratio_thresholds(1:3) > 0)
	error 'The three first logratio thresholds should be negative.';
end
if any(abs_thresholds(1:3) > 0)
	error 'The three first absolute thresholds should be negative.';
end

valid = all(~isnan(test.mean), 2) & all(~isnan(ref.mean), 2);
test = filter_rows(test, valid);
ref = filter_rows(ref, valid);

log_test = log2(test.mean + 1);
log_ref = log2(ref.mean + 1);

ref_min = quantile(log_ref, 1-ref_percentile, 2);
ref_max = quantile(log_ref, ref_percentile, 2);
%ref_min = min(log_ref, [], 2);
%ref_max = max(log_ref, [], 2);

aberrated = zeros(size(test.mean));
for k = 1:3
	aberrated = aberrated - ...
		(log_test - repmat(ref_min, 1, size(log_test, 2)) < ...
		logratio_thresholds(k) & ...
		2.^log_test - repmat(2.^ref_min, 1, size(log_test, 2)) < ...
		abs_thresholds(k));
end
for k = 4:6
	aberrated = aberrated + ...
		(log_test - repmat(ref_max, 1, size(log_test, 2)) > ...
		logratio_thresholds(k) & ...
		2.^log_test - repmat(2.^ref_max, 1, size(log_test, 2)) > ...
		abs_thresholds(k));
end	

