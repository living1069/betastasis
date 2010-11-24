function probes = read_probes_agilent(sample_file)

fid = fopen(sample_file);
while 1
	line = fgetl(fid);
	if line == -1, break, end
	if length(line) < 8, continue, end
	
	if strcmp(line(1:8), 'FEPARAMS')
		[rows, cols] = parse_uarray_size(fid, line);
	elseif strcmp(line(1:8), 'FEATURES')
		[x, y, seq] = parse_agilent_features(fid, line, rows, cols);
	end
end

fclose(fid);

keep = true(length(seq), 1);
for k = 1:length(seq)
	if length(seq{k}) < 30
		keep(k) = false;
	end
end

probes = struct;
probes.XPos = x(keep);
probes.YPos = y(keep);
probes.Sequence = seq(keep);
return;




function [rows, cols] = parse_uarray_size(fid, header)

rows_column = 0;
cols_column = 0;

% Figure out which columns we should read.
column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

for k = 1:length(column_names)
	if strcmp(column_names{k}, 'Grid_NumRows'), rows_column = k; end
	if strcmp(column_names{k}, 'Grid_NumCols'), cols_column = k; end
end

if cols_column == 0 || rows_column == 0
	error 'Invalid file format, vital columns missing.';
end

% Construct a format string for parsing.
parse_format = '';
for k = 1:length(column_names)
	if k == rows_column || k == cols_column
		parse_format = [parse_format '%d'];
	else
		parse_format = [parse_format '%*s'];
	end
end

line = fgetl(fid);
if line == -1, error 'Invalid microarray attributes section.'; end
	
data = textscan(line, parse_format, 'Delimiter', '\t', 'BufSize', 16384);
rows = data{1};
cols = data{2};

return;



function [x, y, seq] = parse_agilent_features(fid, header, rows, cols)

col_column = 0;
row_column = 0;
seq_column = 0;
	
% Figure out which columns we should read.
column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

for k = 1:length(column_names)
	if strcmp(column_names{k}, 'Col'), col_column = k; end
	if strcmp(column_names{k}, 'Row'), row_column = k; end
	if strcmp(column_names{k}, 'Sequence'), seq_column = k; end
end

if col_column == 0 || row_column == 0 || seq_column == 0
	error 'Invalid file format, critical columns missing.';
end


%[~, idx] = sort([row_column, col_column, g_intensity_column, ...
%	r_intensity_column]);
%if sum(idx ~= 1:4) > 0
%	error 'ERROR: Columns are in an unsupported order.';
%end

%fprintf(1, 'Detected feature section with %d columns.\n', length(column_names));
%fprintf(1, 'Probe X and Y positions stored in columns %d and %d.\n', ...
%	col_column, row_column);
%fprintf(1, 'Probe green intensity medians stored in column %d.\n', ...
%	g_intensity_column);
%fprintf(1, 'Probe red intensity medians stored in column %d.\n', ...
%	r_intensity_column);

% Construct a format string for parsing.
parse_format = '';
for k = 1:length(column_names)
	if k == col_column || k == row_column
		parse_format = [parse_format '%d'];
		continue;
	end
	
	if k == seq_column
		parse_format = [parse_format '%s'];
		continue;
	end

	parse_format = [parse_format '%*s'];
end

data = textscan(fid, parse_format, 'Delimiter', '\t', 'BufSize', 16384);
y = data{1};
x = data{2};
seq = data{3};

return;

