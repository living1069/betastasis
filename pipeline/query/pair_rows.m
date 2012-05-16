
% Author: Matti Annala <matti.annala@tut.fi>

function [A, B] = pair_rows(A, B, key)

eval(['keys_a = A.rows.' key ';']);
eval(['keys_b = B.rows.' key ';']);

if length(unique(keys_a)) ~= length(keys_a)
	error 'Row keys in first dataset are not unique.';
end
if length(unique(keys_b)) ~= length(keys_b)
	error 'Row keys in second dataset are not unique.';
end

[common, ia, ib] = intersect(keys_a, keys_b);
A = filter_rows(A, ia);
B = filter_rows(B, ib);

