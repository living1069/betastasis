
% Author: Matti Annala <matti.annala@tut.fi>

function s = nanfilter_rows(s, selected)

fields = fieldnames(s);
fields = fields(~strcmp(fields, 'meta'));    % Leave this field alone.
for k = 1:length(fields)
	f = getfield(s, fields{k});
	
	if isstruct(f)
		s = setfield(s, fields{k}, nanfilter_rows(f, selected));
	elseif ischar(f)
		continue;
	else
		nans = isnan(selected);
		if isnumeric(f)
			v = nan(length(selected), size(f, 2));
		elseif iscellstr(f)
			v = repmat({''}, length(selected), size(f, 2));
		elseif iscell(f)
			v = cell(length(selected), size(f, 2));
		end
		v(~nans, :) = f(selected(~nans), :);
		s = setfield(s, fields{k}, v);
	end
end


