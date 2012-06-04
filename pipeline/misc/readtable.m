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

[data, headers] = readtable_tlm(uncompressed, varargin{:});






function [data, headers] = readtable_tlm(file, varargin)

has_header = true;
num_fields = {};
ignore_fields = {};
comment_regex = '';
num_lines = Inf;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'header')
		if varargin{k+1} > 1
			error 'Field "HeaderLines" is deprecated. Use "Comment" instead.';
		end
		has_header = varargin{k+1};
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
		comment_regex = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, '(num.*lines|lines.*to.*read)')
		num_lines = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end


if ischar(file)
	fid = fopen(file);
else
	fid = file;
end

% If the user has specified a regexp for commented lines, we look for the
% header line after skipping all comment lines. After this code has run,
% the variable "line" will hold the first header/data line. We can't use fseek()
% to simplify the code, because we may need to read from a UNIX pipe.
line = fgetl(fid);
if ~isempty(comment_regex)
	while ischar(line) && rx(line, comment_regex), line = fgetl(fid); end
end

cols = textscan(line, '%s', 'Delimiter', '\t'); cols = cols{1};
if has_header
	headers = cols;
else
	headers = repmat({''}, 1, length(cols));
end

numeric = false(length(headers), 1);
for r = 1:length(num_fields)
	numeric = numeric | rx(headers, num_fields{r});
end

ignore = false(length(headers), 1);
if iscellstr(ignore_fields)
	for r = 1:length(ignore_fields)
		ignore = ignore | rx(headers, ignore_fields{r});
	end
elseif isnumeric(ignore_fields)
	ignore_fields(ignore_fields > length(ignore)) = [];
	ignore(ignore_fields) = true;
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

if num_lines < Inf
	data = textscan(fid, format_str, num_lines, ...
		'Delimiter', '\t', 'BufSize', 16384, 'ReturnOnError', false);
else
	data = textscan(fid, format_str, ...
		'Delimiter', '\t', 'BufSize', 16384, 'ReturnOnError', false);
end

headers = headers(~ignore);
numeric = numeric(~ignore);

% If there was no header, then the first row of data is still stored in the
% variable "cols", and we need to add it to the textscan() parsed data.
if ~has_header
	for k = 1:length(headers)
		if numeric(k)
			data{k} = [str2double(cols{k}); data{k}];
		else
			data{k} = [cols(k); data{k}];
		end
	end
end

if ischar(file), fclose(fid); end

