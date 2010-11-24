
% This function determines specific hybridization affinities for microarray
% probesets by comparing the microarray expression profiles against
% sequencing-based expression profiles. The assumption is that sequencing-based
% profiles are more accurate and have less sequence specific bias.
%
% The two expression profiles must have the same amount of rows.
%
% Inputs:
%     uarray_expr - Microarray-based expression profile.
%     seq_expr - Sequencing-based expression profile.
% 
% Outputs:
%     affinities - Column vector of specific hybridization affinities.
%
% Author: Matti Annala <matti.annala@tut.fi>

function affinities = learn_probeset_affinities(uarray_expr, seq_expr)

if any(size(uarray_expr) ~= size(seq_expr))
	error(['Expression matrices must have the same dimensions when ' ...
	       'calculating the specific hybridization affinities.']);
end

% This is very simple at the moment, but can be extended later.
affinities = median(uarray_expr ./ seq_expr, 2);

