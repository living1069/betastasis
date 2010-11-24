function correlations = robust_pairwise_corr(a, b, corr_type)

type = 'Pearson';
if nargin == 3
	type = corr_type;
end

if size(a, 2) ~= size(b, 2)
	error(['For the calculation of pairwise correlation, matrices a and b ' ...
	       'must have an equal amount of columns.']);
end

correlations = zeros(1, size(a, 2));
for k = 1:size(a, 2)
	d = robust_corr([a(:, k), b(:, k)], type);
	correlations(k) = d(1, 2);
end

