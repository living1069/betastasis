function [samples, array_files, channels] = read_tcga_sdrf(filepath)

fid = fopen(filepath);
header = fgetl(fid);

col_sample = 0;
col_array_file = 0;
col_label = 0;

column_names = regexp(header, '(.+?)(\t|$)', 'tokens');
for k = 1:length(column_names)
	column_name = column_names{k};
	column_name = column_name{1};
	if strcmpi(column_name, 'Extract Name'), col_sample = k; end
	if strcmpi(column_name, 'Array Data File'), col_array_file = k; end
	if strcmpi(column_name, 'Label'), col_label = k; end
end

% Construct a format string for parsing.
parse_format = '';
for k = 1:length(column_names)
	parse_format = [parse_format '%s'];
end

data = textscan(fid, parse_format, 'Delimiter', '\t');
fclose(fid);

samples = data{col_sample};
array_files = data{col_array_file};

empty = strcmp('', samples);
samples = samples(~empty);
array_files = array_files(~empty);

% If the SDRF file has a "Label" column with Cy5 and Cy3 labels, we hold
% on to the channel information.
if col_label && any(strcmpi('Cy5', data{col_label})) && ...
	any(strcmpi('Cy3', data{col_label}))
	channels = data{col_label};
	channels = channels(~empty);
	channels = regexprep(channels, 'cy(\d)', 'Cy$1', 'ignorecase');
else
	channels = {};
end

%[~, uniq] = unique(samples);
%samples = samples(uniq);
%array_files = array_files(uniq);

