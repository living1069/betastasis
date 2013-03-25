
function data = normalize_quantiles(data)

S = size(data, 2);

if any(any(isnan(data))) || any(any(data == Inf))
	error('Data contains NaN or Inf values.');
end

%sketch_cols = 1:floor((size(samples.mean, 2)-1)/S):size(samples.mean, 2);
sketch_mean = zeros(size(data, 1), 1);

progress = Progress;
for s = 1:S
	[sorted, data(:, s)] = sort(data(:, s), 'descend');
	sketch_mean = sketch_mean + sorted / S;
	progress.update(0.7 * s/S);
end
for s = 1:S
	data(data(:, s), s) = sketch_mean;
	progress.update(0.7 + 0.3 * s/S);
end


