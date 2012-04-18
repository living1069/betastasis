% RX    Run regular expression on cell array of strings, return logical.
%
%    MATCHES = RX(X, REGEX)

% Author: Matti Annala <matti.annala@tut.fi>

function m = rx(x, regex)

if ischar(x), x = { x }; end

m = false(size(x));
for k = 1:numel(x)
	if regexpi(x{k}, regex), m(k) = true; end
end

