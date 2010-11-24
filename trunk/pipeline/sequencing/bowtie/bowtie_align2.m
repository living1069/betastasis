function results = bowtie_align(reads, index, varargin)

global pipeline_config;

max_mismatches = 0;
read_len_filter = [];
trim_len = NaN;
report_alignments = [];
allow_alignments = [];
max_threads = pipeline_config.MaxThreads;
columns = 'read,strand,target,offset,sequence,quality';
output_file = '';

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MaxMismatches')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Trim')
		trim_len = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'LengthFilter')
		read_len_filter = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MaxThreads')
		max_threads = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'ReportAlignments')
		report_alignments = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'AllowAlignments')
		allow_alignments = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Columns')
		columns = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'OutputFile')
		output_file = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

[flags, index_suffix] = bowtie_flags_for_reads(reads);

if ischar(index)
	% If the index name is specified as a relative path, prefix it with the
	% pipeline directory that holds all Bowtie indices.
	if index(1) ~= '/'
		index = bowtie_index(index);
	end
	index_name = [index index_suffix];
else
	[color, ~] = seq_read_type(reads);
	fasta_tmp = ptemp;
	index_tmp = ptemp;

	write_seq_fasta(index, fasta_tmp);

	color_option = '';
	index_suffix = '';
	if color
		color_option = '-C';
		index_suffix = '_colorspace';
	end

	index_name = [index_tmp index_suffix];

	[status, ~] = unix(sprintf('%s/tools/bowtie/bowtie-build %s %s %s', ...
		ppath, color_option, fasta_tmp, index_name));
	safe_delete(fasta_tmp);
	if status ~= 0, error 'Bowtie index construction failed.'; end
end

flags = [flags ' -B1'];
flags = sprintf('%s -v%d', flags, max_mismatches);
flags = sprintf('%s -p%d', flags, max_threads);

bcols = [];
tokens = regexp(columns, '(\w+?)(,|\s|$)', 'tokens');
for k = 1:length(tokens)
	token = tokens{k}; token = token{1};
	if regexpi(token, 'read|readid|read_id|id'), bcols(end+1) = 1;
	elseif regexpi(token, 'strand'), bcols(end+1) = 2;
	elseif regexpi(token, 'ref|target'), bcols(end+1) = 3;
	elseif regexpi(token, 'offset|pos|position'), bcols(end+1) = 4;
	elseif regexpi(token, 'sequence|seq'), bcols(end+1) = 5;
	elseif regexpi(token, 'quality|qual'), bcols(end+1) = 6;
	elseif regexpi(token, 'mismatches'), bcols(end+1) = 8;
	else error('Invalid column identifier "%s".', token);
	end
end

suppress_cols = sprintf('%d,', setdiff(1:8, bcols));
suppress_cols = suppress_cols(1:end-1);
flags = sprintf('%s --suppress %s', flags, suppress_cols);

if ~isempty(report_alignments) && ~isempty(allow_alignments)
	fprintf(1, ['WARNING: Specified AllowAlignments and ReportAlignments. ' ...
		'Using only the AllowAlignments option.\n']);
	report_alignments = [];
end

if ~isempty(report_alignments)
	flags = sprintf('%s -k%d', flags, report_alignments);
end

if ~isempty(allow_alignments)
	flags = sprintf('%s -m%d', flags, allow_alignments);
end

if ~isempty(read_len_filter)
	filtered_tmp = filter_reads_by_len(reads, read_len_filter);
else
	filtered_tmp = reads;
end

if trim_len > 0
	trimmed_reads_tmp = ptemp;
	
	[color, quality] = seq_read_type(reads);
	
	fprintf(1, 'Trimming reads to a length of %d bp...\n', trim_len);
	[status, ~] = unix(sprintf( ...
		'%s/sources/sequencing/trim_reads.py %s %d > %s', ...
		ppath, filtered_tmp, trim_len + 3 * color, trimmed_reads_tmp));
else
	trimmed_reads_tmp = filtered_tmp;
end

fprintf(1, 'Using flags "%s" when invoking Bowtie.\n', flags);

if isempty(output_file)
	alignments_file = ptemp;
else
	alignments_file = output_file;
end

[status, out] = unix(sprintf('%s/tools/bowtie/bowtie %s %s %s > %s', ...
	ppath, flags, index_name, trimmed_reads_tmp, alignments_file));
if status ~= 0
	fprintf(1, '%s', out);
	error 'Bowtie read alignment failed.';
end

results = struct;

out = strread(out, '%s', 'delimiter', '\n');
for k = 1:length(out)
	match = regexp(out{k}, '^# reads processed: (\d+)', 'tokens');
	if length(match) == 1
		tmp = match{1};
		results.TotalReads = str2double(tmp{1});
	end
end
	
if ~isempty(read_len_filter)
	safe_delete(filtered_tmp);
end
	
if trim_len > 0
	safe_delete(trimmed_reads_tmp);
end

% If the user wants to keep the output file, then we assume that he is not
% interested in the parsed file contents.
if ~isempty(output_file)
	return;
end

parse_format = '';
file_col = 0;
if ismember(1, bcols)
	parse_format = [parse_format '%s'];
end
if ismember(2, bcols)
	parse_format = [parse_format '%s'];
end
if ismember(3, bcols)
	parse_format = [parse_format '%s'];
end
if ismember(4, bcols)
	parse_format = [parse_format '%d'];
end
if ismember(5, bcols)
	parse_format = [parse_format '%s'];
end
if ismember(6, bcols)
	parse_format = [parse_format '%s'];
end
if ismember(8, bcols)
	parse_format = [parse_format '%s'];
end

fid = fopen(alignments_file);
data = textscan(fid, parse_format);
fclose(fid);

delete(alignments_file);

file_col = 0;
if ismember(1, bcols)
	file_col = file_col + 1; results.ReadID = data{file_col};
end
if ismember(2, bcols)
	file_col = file_col + 1;
	strands = data{file_col};
	results.Strand = repmat(' ', length(strands), 1);
	for k = 1:length(strands)
		results.Strand(k) = strands{k};
	end
end
if ismember(3, bcols)
	file_col = file_col + 1; results.Target = data{file_col};
end
if ismember(4, bcols)
	file_col = file_col + 1; results.Offset = data{file_col};
end
if ismember(5, bcols)
	file_col = file_col + 1; results.Sequence = data{file_col};
end
if ismember(6, bcols)
	file_col = file_col + 1; results.Quality = data{file_col};
end
if ismember(8, bcols)
	file_col = file_col + 1; results.Mismatches = data{file_col};
end












function [flags, index_suffix] = bowtie_flags_for_reads(reads)

[color, quality] = seq_read_type(reads);

flags = '';
index_suffix = '';

if color == 1
	index_suffix = '_colorspace';
end

if color == 1 && quality == 1
	flags = '-C --solexa-quals --integer-quals';
elseif color == 1 && quality == 0
	flags = '-C -f';
elseif color == 0 && quality == 1
	flags = '';
elseif color == 0 && quality == 0
	flags = '-f';
end

