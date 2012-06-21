
% Author: Matti Annala <matti.annala@tut.fi>

function m = rx_capture(x, regex)

if ischar(x), x = { x }; end
tokens = regexpi(x, regex, 'tokens');
if any(cellfun(@isempty, tokens))
	error 'Regexp capture failed for some items.';
end
tokens = vertcat(tokens{:});
m = vertcat(tokens{:});

