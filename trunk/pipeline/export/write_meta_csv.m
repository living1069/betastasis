
% Author: Matti Annala <matti.annala@tut.fi>

function [] = write_meta_csv(meta, filepath)

sheet = {};
fields = fieldnames(meta);

while ~isempty(fields)
	fname = fields{1};
	fields = fields(2:end);
	
	f = eval(sprintf('meta.%s', fname));
	
	if iscellstr(f)
		sheet{1, end+1} = fname;
		sheet(2:length(f)+1, end) = f;
		
	elseif isnumeric(f) || islogical(f)
		sheet{1, end+1} = fname;
		strs = arrayfun(@num2str, f, 'UniformOutput', false);
		nan_strs = strcmp(strs, 'NaN');
		strs(nan_strs) = repmat({''}, 1, sum(nan_strs));
		sheet(2:length(f)+1, end) = strs;
		
	elseif isstruct(f)
		fields = [fields; strcat([fname '.'], fieldnames(f))];
		
	else
		error('Unsupported metadata field "%s".', fields{f});
	end
end

fid = fopen(filepath, 'W');
for y = 1:size(sheet, 1)
	fprintf(fid, '%s\t', sheet{y, 1:end-1});
	fprintf(fid, '%s\n', sheet{y, end});
end
fclose(fid);



