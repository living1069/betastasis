
% Author: Matti Annala <matti.annala@tut.fi>

function filtered = filter(ds, selected)

if islogical(selected) || ...
	(isnumeric(selected) && ~isempty(selected) && min(selected) > 0) || ...
	isempty(selected)
else
	error 'filter() must be called with logical or index vector.';
end

filtered = filter_struct(ds, selected);





function s = filter_struct(s, selected)

fields = fieldnames(s);
fields = fields(~strcmp(fields, 'rows'));    % Leave this field alone.
for k = 1:length(fields)
	f = getfield(s, fields{k});
	if isstruct(f)
		s = setfield(s, fields{k}, filter_struct(f, selected));
	elseif ischar(f)
		continue;
	else
		s = setfield(s, fields{k}, f(:, selected));
	end
end

