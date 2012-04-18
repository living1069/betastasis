function [data, headers] = readtable(file, varargin)

if rx(file, '.*\.gz$')
	tmp = temporary('readtable_compressed');
	
	uncompressed = [tmp 'uncompress_pipe'];
	unix(sprintf('mkfifo %s', uncompressed));
	unix(sprintf('gunzip -c %s > %s &', file, uncompressed));
else
	uncompressed = file;
end

if rx(file, '.*\.xls(.gz)?$')
	[data, headers] = readtable_xls(uncompressed, varargin{:});
else
	[data, headers] = readtable_tlm(uncompressed, varargin{:});
end






function [data, headers] = readtable_tlm(file, varargin)

header_lines = 1;
num_fields = {};

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'header.*line')
		header_lines = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Numeric')
		num_fields = varargin{k+1};
		if ischar(num_fields), num_fields = { num_fields }; end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end


if ischar(file)
	fid = fopen(file);
else
	fid = file;
end

if header_lines > 0
	for k = 1:header_lines-1, fgetl(fid); end
	header = fgetl(fid);
	cols = textscan(header, '%s', 'Delimiter', '\t');
	headers = cols{1};
elseif header_lines == 0
	line = fgetl(fid);
	cols = textscan(line, '%s', 'Delimiter', '\t');
	headers = repmat({''}, 1, length(cols{1}));
	frewind(fid);
end

format_str = '';
numeric = false(length(headers), 1);
for r = 1:length(num_fields)
	numeric = numeric | rx(headers, num_fields{r});
end

for k = 1:length(headers)
	if numeric(k)
		format_str = [format_str '%f'];
	else
		format_str = [format_str '%s'];
	end
end

format_str = [format_str '%*s'];

data = textscan(fid, format_str, 'Delimiter', '\t', 'BufSize', 16384, ...
	'ReturnOnError', false);
	
if ischar(file), fclose(fid); end

	
	
	
	
	
	
function [data, headers] = readtable_xls(file, varargin)

[~, ~, xls] = xlsread(file)

