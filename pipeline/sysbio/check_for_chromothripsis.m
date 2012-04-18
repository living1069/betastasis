
% Author: Matti Annala <matti.annala@tut.fi>

function [] = check_for_chromothripsis(segments)

num_segs = nan(size(segments.chromosome));
for k = 1:numel(segments.chromosome)
	num_segs(k) = length(segments.chromosome{k}.logratio);
end

num_segs

