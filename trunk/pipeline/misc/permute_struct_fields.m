function permuted = permute_struct_fields(s, indices)

permuted = struct;
unknown = isnan(indices);
indices(unknown) = 1;

fields = fieldnames(s);
for k = 1:length(fields)
	f = getfield(s, fields{k});
	if isnumeric(f) || iscell(f)
		nf = f(indices, :);
		if isnumeric(f), nf(unknown) = NaN;
		elseif iscellstr(f), nf(unknown) = repmat({'-'}, sum(unknown), 1);
		elseif iscell(f), nf(unknown) = repmat({[]}, sum(unknown), 1); 
		end
		permuted = setfield(permuted, fields{k}, nf);
	else
		error(['permute_struct_fields() does not know how to handle ' ...
		       'field ''%s''.\n'], fields{k});
	end
end

