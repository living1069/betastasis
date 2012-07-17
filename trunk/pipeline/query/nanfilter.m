
% Author: Matti Annala <matti.annala@tut.fi>

function filtered = nanfilter(ds, selected)

if ~isnumeric(selected)
	error 'nanfilter() must be called with an index vector.';
end

filtered = filter_columns(ds, selected);





function s = filter_columns(s, selected)

fields = fieldnames(s);
fields = fields(~strcmp(fields, 'rows'));    % Leave this field alone.
for k = 1:length(fields)
	f = getfield(s, fields{k});
	
	if isstruct(f)
		s = setfield(s, fields{k}, filter_columns(f, selected));
	elseif ischar(f)
		continue;
	else
		nans = isnan(selected);
		if isnumeric(f)
			v = nan(size(f, 1), length(selected));
		elseif iscellstr(f)
			v = repmat({''}, size(f, 1), length(selected));
		elseif iscell(f)
			v = cell(size(f, 1), length(selected));
		end
		v(:, ~nans) = f(:, selected(~nans));
		s = setfield(s, fields{k}, v);
	end
end

