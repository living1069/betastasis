% RX    Run regular expression on cell array of strings, return logical.
%
%    MATCHES = RX(X, REGEX)

% Author: Matti Annala <matti.annala@tut.fi>

function m = rx(x, regex)

if ischar(x), x = { x }; end
m = ~cellfun(@isempty, regexpi(x, regex));

