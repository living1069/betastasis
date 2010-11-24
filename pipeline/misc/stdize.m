function x = stdize(x)
for k = 1:size(x, 2)
	x(:, k) = x(:, k) - mean(x(:, k));
	x(:, k) = x(:, k) / std(x(:, k));
end

