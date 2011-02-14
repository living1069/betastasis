function data = read_uarray_sample_agilent(sample_file, probes, channels)

N = length(probes.Sequence);

if nargin < 3, channels = {}; end

fid = fopen(sample_file);

% Sanity check, see if the file even is an Agilent array file.
line = fgetl(fid);
if ~strcmp(line(1:4), 'TYPE')
	error 'The given file does not look like an Agilent sample file.';
end

while 1
	line = fgetl(fid);
	if line == -1, break, end
	if length(line) < 8, continue, end
	
	if strcmp(line(1:8), 'FEPARAMS')
		[rows, cols] = parse_uarray_size(fid, line);
	elseif strcmp(line(1:8), 'FEATURES')
		cdata = parse_agilent_features(fid, line, rows, cols, channels);
	end
end

if ~exist('cdata', 'var')
	error 'Could not find a FEATURES section in the file.';
end

fclose(fid);

data = zeros(N, length(channels));
for k = 1:length(cdata)
	raw = cdata{k};
	for p = 1:N
		data(p, k) = raw(probes.YPos(p), probes.XPos(p));
	end
end

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

%fprintf(1, 'Microarray has %d rows and %d columns.\n', rows, cols);
return;



function channel_data = parse_agilent_features(fid, header, rows, cols, ...
	channels)

col_column = 0;
row_column = 0;
channel_columns = [];
	
% Figure out which columns we should read.
column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

for k = 1:length(column_names)
	if strcmp(column_names{k}, 'Col'), col_column = k; end
	if strcmp(column_names{k}, 'Row'), row_column = k; end
end

if col_column == 0 || row_column == 0
	error 'Invalid file format, probe position columns missing.';
end

if isempty(channels)
	if any(strcmpi('gProcessedSignal', column_names))
		channels = { 'gProcessedSignal' };
		if any(strcmpi('rProcessedSignal', column_names))
			channels{2} = 'rProcessedSignal';
		end
	else
		error(['User did not specify channel names and cannot find any ' ...
		       'known channels.']);
	end
end

if length(channels) > 0
	for k = 1:length(column_names)
		for c = 1:length(channels)
			if strcmpi(column_names{k}, channels{c})
				channel_columns(c) = k;
			end
		end
	end
end

for k = 1:length(channels)
	if channel_columns(k) == 0
		error('Could not find requested column "%s".', channels{k});
	end
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
	
	found = false;
	for c = 1:length(channels)
		if k == channel_columns(c)
			parse_format = [parse_format '%f'];
			found = true;
			break;
		end
	end
	
	if ~found
		parse_format = [parse_format '%*s'];
	end
end

data = textscan(fid, parse_format, 'Delimiter', '\t', 'BufSize', 16384);
probe_y = data{1};
probe_x = data{2};

[~, channel_order] = sort(channel_columns);

channel_data = cell(length(channels), 1);

for c = 1:length(channels)
	channel = channel_order(c);
	raw = zeros(rows, cols);
	channel_raw = data{2+c};
	
	for k = 1:length(channel_raw)
		raw(probe_y(k), probe_x(k)) = channel_raw(k);
	end
	
	channel_data{c} = raw;
end

return;
