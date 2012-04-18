function p = ttest_mean_std(a_mean, a_std, a_n, b_mean, b_std, b_n)

t = (a_mean - b_mean) / sqrt(a_std^2 / a_n + b_std^2 / b_n);

a_sn = a_std^2 / a_n;
b_sn = b_std^2 / b_n;
df = (a_sn + b_sn)^2 / ((a_sn^2 / (a_n-1)) + (b_sn^2) / (b_n - 1));

% Two-tailed.
p = tcdf(t, df);
if p > 0.5, p = 1 - p; end
p = p * 2;

