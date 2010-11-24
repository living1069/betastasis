function s = filter_struct(s, selected)

if ~islogical(selected) && ~(isnumeric(selected) && min(selected) > 0)
	error 'filter_struct() called with an invalid second argument.';
end

fields = fieldnames(s);
for k = 1:length(fields)
	f = getfield(s, fields{k});
	if isstruct(f)
		s = setfield(s, fields{k}, filter_struct(f, selected));
	elseif ischar(f)
		continue;
	else
		s = setfield(s, fields{k}, f(selected, :));
	%else
	%	error('filter_struct() found a field "%s" it cannot handle.\n', ...
	%		fields{k});
	end
end

