
% ASSOCIATE   Look for co-occurring or mutually exclusive events 
%
%    [PVAL, POSITIVE] = ASSOCIATE(X, Y) uses the hypergeometric distribution
%    to calculate p-values for the co-occurrence or mutual exclusivity of 
%    events described in the logical column vector Y and the events described
%    in each column of logical matrix X.

% Author: Matti Annala <matti.annala@tut.fi>

function [pval, positive] = associate(X, y)

M = length(y);
K = sum(X, 1);
N = sum(y);

x = sum(X & repmat(y, 1, size(X, 2)), 1);
p_exclusive = hygecdf(x, M, K, N);

x = sum(X & repmat(y, 1, size(X, 2)), 1) - 1;
p_mutual = 1 - hygecdf(x, M, K, N);

pval = nan(1, size(X, 2));
positive = p_mutual < p_exclusive;
pval(positive) = p_mutual(positive);
pval(~positive) = p_exclusive(~positive);

