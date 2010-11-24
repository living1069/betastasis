function correlation = robust_corr(a, corr_type)

type = 'Pearson';
if nargin == 2
	type = corr_type;
end

b = sum(a, 2);
bad_nan = isnan(b);
bad_ninf = (b == -Inf);
bad_inf = (b == Inf);
bad = bad_nan | bad_inf;

% Replace all negative infinities with the smallest finite value in the data.
ninf = (a == -Inf);
a(ninf) = Inf;
a(ninf) = min(min(a));

fprintf(1, ['Robust correlation calculated for %.1f%% of all rows ' ...
            '(NaN %.1f%%, -Inf %.1f%%, Inf %.1f%%)\n'], ...
	sum(~bad) / length(b) * 100, sum(bad_nan) / length(b) * 100, ...
	sum(bad_ninf) / length(b) * 100, sum(bad_inf) / length(b) * 100);

correlation = corr(a(~bad, :), 'type', type);

