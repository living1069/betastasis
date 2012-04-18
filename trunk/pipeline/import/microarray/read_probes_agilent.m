function probes = read_probes_agilent(sample_file, sequence_map)

fid = fopen(sample_file);
while 1
	line = fgetl(fid);
	if line == -1, break, end
	if length(line) < 8, continue, end
	
	if strcmp(line(1:8), 'FEPARAMS')
		[rows, cols] = parse_uarray_size(fid, line);
	elseif strcmp(line(1:8), 'FEATURES')
		[x, y, seq, probe_id] = parse_agilent_features(fid, line, rows, cols);
	end
end

fclose(fid);


% If the file did not contain any sequence information, we require that the
% user must specify a mapping from probe IDs to probe sequences.
if isempty(seq)
	valid = sequence_map.isKey(probe_id);
	if any(~valid)
		fprintf(1, 'Discarded %d unknown probe IDs. Some examples:\n', ...
			sum(~valid));
		examples = find(~valid);
		if length(examples) > 5, examples = examples(1:5); end
		fprintf(1, '%s\n', probe_id{examples});
	end

	seq = repmat({''}, length(valid), 1);
	seq(valid) = sequence_map.values(probe_id(valid));
end

% Check that the probe sequences are reasonable.
whos valid seq
valid = rx(seq, '[ACGT]+');
	
if any(~valid)
	fprintf(1, 'Discarded %d probes in total.\n', sum(~valid));
end

probes = struct;
probes.xpos = x(valid);
probes.ypos = y(valid);
probes.sequence = seq(valid);







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








function [x, y, seq, probe_id] = parse_agilent_features(fid, header, rows, cols)

col_column = 0;
row_column = 0;
seq_column = 0;
probe_column = 0;
	
% Figure out which columns we should read.
column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

for k = 1:length(column_names)
	if strcmp(column_names{k}, 'Col'), col_column = k; end
	if strcmp(column_names{k}, 'Row'), row_column = k; end
	if strcmp(column_names{k}, 'Sequence'), seq_column = k; end
	if rx(column_names{k}, 'probename'), probe_column = k; end
end

if col_column == 0 || row_column == 0 || (seq_column == 0 && probe_column == 0)
	error 'Invalid file format, critical columns missing.';
end

if seq_column > 0 && probe_column > 0
	probe_column = 0;
end

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
	
	if k == seq_column || k == probe_column
		parse_format = [parse_format '%s'];
		continue;
	end
	
	parse_format = [parse_format '%*s'];
end

data = textscan(fid, parse_format, 'Delimiter', '\t', 'BufSize', 16384);
y = data{1};
x = data{2};

probe_id = {};
seq = {};
if probe_column > 0
	probe_id = data{3};
else
	seq = data{3};
end

