
% Author: Matti Annala <matti.annala@tut.fi>

function s = filter_rows(s, selected)

fields = fieldnames(s);
fields = fields(~strcmp(fields, 'meta'));    % Leave this field alone.
for k = 1:length(fields)
	f = getfield(s, fields{k});
	if isstruct(f)
		s = setfield(s, fields{k}, filter_rows(f, selected));
	elseif ischar(f)
		continue;
	else
		s = setfield(s, fields{k}, f(selected, :));
	end
end



