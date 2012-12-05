
% ASSOCIATE   Look for co-occurring or mutually exclusive events 
%
%    [PVAL, POSITIVE] = ASSOCIATE(X, Y) uses the hypergeometric distribution
%    to calculate p-values for the co-occurrence or mutual exclusivity of 
%    events described in the logical column vector Y and the events described
%    in each column of logical matrix X.
%    
%    NaN values are treated as missing values.

% Author: Matti Annala <matti.annala@tut.fi>

function [pval, positive] = associate(X, y)

M = sum(~isnan(X), 1);   % How many patients in total?
K = nansum(X, 1);        % How many mutation positive patients?

% How many tumor patients?
N = nansum(~isnan(X) & repmat(y, 1, size(X, 2)), 1);

% How many mutation positive tumor patients?
x = nansum((X > 0) & repmat(y, 1, size(X, 2)), 1);  

p_exclusive = hygecdf(x, M, K, N);

x = nansum((X > 0) & repmat(y, 1, size(X, 2)), 1) - 1;
p_mutual = 1 - hygecdf(x, M, K, N);

pval = nan(1, size(X, 2));
positive = p_mutual < p_exclusive;
pval(positive) = p_mutual(positive);
pval(~positive) = p_exclusive(~positive);

