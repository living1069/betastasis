function [] = robust_scatter(a, b, varargin)

bad_nan = isnan(a) | isnan(b);
bad_ninf = (a == -Inf) | (b == -Inf);
bad_inf = (a == Inf) | (b == Inf);
bad = bad_nan | bad_ninf | bad_inf;

fprintf(1, ['Robust scatter plot shown for %.1f%% of all pairs ' ...
            '(NaN %.1f%%, -Inf %.1f%%, Inf %.1f%%)\n'], ...
	sum(~bad) / length(b) * 100, sum(bad_nan) / length(b) * 100, ...
	sum(bad_ninf) / length(b) * 100, sum(bad_inf) / length(b) * 100);

scatter(a(~bad), b(~bad), varargin{:});

