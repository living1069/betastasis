
% Author: Matti Annala <matti.annala@tut.fi>

function [pval, positive] = associate(X, y)

M = length(y);
K = sum(X, 1);
N = sum(y);

x = sum(X & repmat(y, 1, size(X, 2)), 1);
p_exclusive = hygecdf(x, M, K, N);

x = sum(X & repmat(y, 1, size(X, 2)), 1) - 1;
x = max(x, 0);
p_mutual = 1 - hygecdf(x, M, K, N);

pval = nan(1, size(X, 2));
positive = p_mutual < p_exclusive;
pval(positive) = p_mutual(positive);
pval(~positive) = p_exclusive(~positive);

