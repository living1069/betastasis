function merged = cat_structs(varargin)

if length(varargin) < 2
	error 'cat_structs() requires at least two input arguments.';
end

merged = struct;

fields = {};
for k = 1:length(varargin)
	fields = union(fields, fieldnames(varargin{k}));
end

rows = 0;
new_rows = 0;

% WARNING: Some tricky code ahead.

for p = 1:length(varargin)
	rows = rows + new_rows;
	new_rows = 0;
	
	missing_fields = {};
	
	for k = 1:length(fields)
		if ~isfield(varargin{p}, fields{k})
			% If a field is missing from the structure, we will handle it later
			% once we know how many rows the structure contains in other fields.
			missing_fields{end + 1} = fields{k};
			
		elseif ~isfield(merged, fields{k})
			% If there is a new field in the structure that is about to be
			% appended, we add the new field to the merged structure with some
			% null data at the beginning.
			f = getfield(varargin{p}, fields{k});
			if isstruct(f) || ischar(f)
				merged = setfield(merged, fields{k}, f);
			elseif iscell(f) && size(f, 2) == 1
				merged = setfield(merged, fields{k}, ...
					cat(1, repmat({'-'}, rows, 1), f));
				new_rows = max(new_rows, size(f, 1));
			elseif (isnumeric(f) || islogical(f)) && size(f, 2) == 1
				merged = setfield(merged, fields{k}, cat(1, nan(rows, 1), f));
				new_rows = max(new_rows, size(f, 1));
			else
				error(['cat_structs() does not know how to handle field ' ...
					   '''%s''.\n'], fields{k});
			end
		
		else
			mf = getfield(merged, fields{k});
			f = getfield(varargin{p}, fields{k});
			
			if isstruct(mf) && isstruct(f)
				eval(['merged.' fields{k} ' = cat_structs(mf, f);']);
			elseif size(mf, 2) == 1
				eval(['merged.' fields{k} ' = cat(1, mf, f);']);
				new_rows = max(new_rows, size(f, 1));
			elseif ischar(mf)
				continue;
			else
				error(['cat_structs() does not know how to handle field ' ...
					   '''%s''.\n'], fields{k});
			end
		end
	end
	
	for k = 1:length(missing_fields)
		field = missing_fields{k};

		% Don't do anything if the field is missing from both the merged
		% structure and the structure that is being appended.
		if ~isfield(merged, field), continue, end

		mf = getfield(merged, field);

		if iscell(mf) && size(mf, 2) == 1
			eval(['merged.' field ' = cat(1, merged.' field ...
				', repmat({''-''}, new_rows, 1));']);
		elseif isnumeric(mf) && size(mf, 2) == 1
			eval(['merged.' field ' = cat(1, merged.' field ...
				', nan(new_rows, 1));']);
		elseif ischar(mf)
			continue;
		else
			error(['cat_structs() does not know how to handle field ' ...
				   '''%s''.\n'], field);
		end
	end
end

