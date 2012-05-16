function [data, headers] = readtable(file, varargin)

if isnumeric(file)
	uncompressed = file;
elseif ischar(file)
	if rx(file, '.*\.gz$')
		tmp = temporary('readtable_compressed');
		
		uncompressed = [tmp 'uncompress_pipe'];
		unix(sprintf('mkfifo %s', uncompressed));
		unix(sprintf('gunzip -c %s > %s &', file, uncompressed));
	else
		uncompressed = file;
	end
end

%if rx(file, '.*\.xls(.gz)?$')
%	[data, headers] = readtable_xls(uncompressed, varargin{:});
%else
	[data, headers] = readtable_tlm(uncompressed, varargin{:});
%end






function [data, headers] = readtable_tlm(file, varargin)

header_lines = 1;
num_fields = {};
ignore_fields = {};
comment_regex = '';

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'header.*line')
		header_lines = varargin{k+1};
		continue;
	end
	
	if rx(varargin{k}, 'numeric')
		num_fields = varargin{k+1};
		if ischar(num_fields), num_fields = { num_fields }; end
		continue;
	end
	
	if rx(varargin{k}, 'ignore')
		ignore_fields = varargin{k+1};
		if ischar(ignore_fields), ignore_fields = { ignore_fields }; end
		continue;
	end
	
	if rx(varargin{k}, 'comment')
		comment_regex = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end


if ischar(file)
	fid = fopen(file);
else
	fid = file;
end

% If the user has specified a regexp for commented lines, we look for the
% header line only after skipping all comment lines.
if ~isempty(comment_regex)
	fpos = ftell(fid);
	comment_lines = 0;
	while 1
		line = fgetl(fid);
		if ~ischar(line) || ~rx(line, comment_regex), break, end
		comment_lines = comment_lines + 1;
	end
	fseek(fid, fpos, -1);
	for k = 1:comment_lines, fgetl(fid); end
end

% If we have header lines, we calculate the number of columns based on them.
if header_lines > 0
	for k = 1:header_lines-1, fgetl(fid); end
	header = fgetl(fid);
	cols = textscan(header, '%s', 'Delimiter', '\t');
	headers = cols{1};
elseif header_lines == 0
	fpos = ftell(fid);
	line = fgetl(fid);
	cols = textscan(line, '%s', 'Delimiter', '\t');
	headers = repmat({''}, 1, length(cols{1}));
	fseek(fid, fpos, -1);
end

numeric = false(length(headers), 1);
for r = 1:length(num_fields)
	numeric = numeric | rx(headers, num_fields{r});
end

ignore = false(length(headers), 1);
for r = 1:length(ignore_fields)
	ignore = ignore | rx(headers, ignore_fields{r});
end

format_str = '';
for k = 1:length(headers)
	if numeric(k)
		format_str = [format_str '%f'];
	elseif ignore(k)
		format_str = [format_str '%*s'];
	else
		format_str = [format_str '%s'];
	end
end

format_str = [format_str '%*s'];

data = textscan(fid, format_str, 'Delimiter', '\t', 'BufSize', 16384, ...
	'ReturnOnError', false);
headers = headers(~ignore);
	
if ischar(file), fclose(fid); end

	
	
	
	
	
	
function [data, headers] = readtable_xls(file, varargin)

[~, ~, xls] = xlsread(file)

