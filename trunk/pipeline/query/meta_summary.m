function [] = meta_summary(meta)

if isfield(meta, 'Meta')
	meta = meta.Meta;
end

print_field(meta, '');





function [] = print_field(field, prefix)

fname = prefix(1:end-1);

if strcmp(fname, 'Resource')
	return;
end

if iscellstr(field)
	fprintf(1, '%s:\n', fname);
	vals = unique(field);
	if length(vals) > 20 
		fprintf(1, '%d different values (not shown).\n', length(vals));
	else
		na_count = 0;
		for k = 1:length(vals)
			if strcmp(vals{k}, '-')
				na_count = sum(strcmp(vals{k}, field));
				continue;
			end
			fprintf(1, '- %s (%d items)\n', vals{k}, ...
				sum(strcmp(vals{k}, field)));
		end
		if na_count
			fprintf(1, '- N/A (%d items)\n', na_count);
		end
	end
	fprintf(1, '\n');
elseif isstruct(field) && numel(field) == 1
	fields = fieldnames(field);
	for k = 1:length(fields)
		eval(['print_field(field.' fields{k} ...
			', [prefix ''' fields{k} '.'']);']);
	end
elseif isnumeric(field) && numel(field) == size(field, 1)
	fprintf(1, '%s:\n', fname);
	fprintf(1, 'Minimum = %.2f\t\tMaximum = %.2f\n', min(field), max(field));
	fprintf(1, 'Mean = %.2f\t\tStdev = %.2f\n',  nanmean(field), nanstd(field));
	fprintf(1, 'N/A count = %d\n', sum(isnan(field)));
	fprintf(1, '\n');
elseif islogical(field) && numel(field) == size(field, 1)
	fprintf(1, '%s:\n', fname);
	fprintf(1, 'True = %d\t\tFalse = %d\t\tN/A = %d\n', sum(field == true), ...
		sum(field == false), sum(isnan(field)));
	fprintf(1, '\n');
end



