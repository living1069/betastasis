function alignments = bowtie_align(reads, index, varargin)

global pipeline_config;

tmp_pool = FilePool;

ignore_pairs = false;
max_mismatches = 0;
report_alignments = [];
allow_alignments = [];
max_threads = pipeline_config.MaxThreads;
columns = 'read,strand,target,offset,sequence,quality';
output_file = '';
unaligned_file = '';

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'IgnorePairs')
		ignore_pairs = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MaxMismatches')
		max_mismatches = varargin{k+1};
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
	
	if strcmpi(varargin{k}, 'UnalignedFile')
		unaligned_file = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

%if isstruct(reads) && isfield(reads.Meta.Sample, 'ID')
%	sample_id = reads.Meta.Sample.ID{1};
%else
%	sample_id = sprintf('sample #%d', r);
%end
fprintf(1, 'Aligning sample using Bowtie:\n');
fprintf(1, '-> Preparing reads for alignment...\n');

paired = false;
color = false;

if isstruct(reads)
	if regexpi(reads.Meta.Sample.SequenceType{1}, 'paired')
		paired = true;
	end
	if regexpi(reads.Meta.Sample.SequenceType{1}, 'color')
		color = true;
	end
end






% Determine general alignment flags.
if color
	flags = '-C -r';
	index_suffix = '_colorspace';
else
	flags = '-r';
	index_suffix = '';
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
	else error('Invalid column identifier "%s".', token);
	end
end

suppress_cols = sprintf('%d,', setdiff(1:8, bcols));
suppress_cols = suppress_cols(1:end-1);
flags = sprintf('%s --suppress %s', flags, suppress_cols);

if report_alignments > allow_alignments
	error 'ReportAlignments should never be higher than AllowAlignments.';
end

if ~isempty(report_alignments)
	flags = sprintf('%s -k%d', flags, report_alignments);
end
if ~isempty(allow_alignments)
	flags = sprintf('%s -m%d', flags, allow_alignments);
end

if ~isempty(unaligned_file)
	flags = sprintf('%s --un %s', flags, unaligned_file);
end





% Determine the type of index that the user wishes to align against.
if ischar(index)
	% If the index name is specified as a relative path, prefix it with the
	% pipeline directory that holds all Bowtie indices.
	if index(1) ~= '/'
		index = bowtie_index(index);
	end
	index_name = [index index_suffix];
else
	fasta_tmp = tmp_pool.temp('index_fasta');
	index_tmp = tmp_pool.temp('index');

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
	if status ~= 0, error 'Bowtie index construction failed.'; end
end




extracted = extract_reads(reads);
read_files = extracted.Raw{1};

if isempty(output_file)
	alignments_file = tmp_pool.temp('alignments');
else
	alignments_file = output_file;
end

al = [];

if paired && ignore_pairs == false
	% The user wishes to perform proper paired end alignment.
	error 'Proper paired end alignment not supported yet.';

elseif paired && ignore_pairs == true
	% The user wishes to align both ends of each pair separately.
	if ~isempty(output_file)
		error 'Output file not supported with paired reads';
	end
	
	fprintf(1, '-> Invoking Bowtie with flags "%s"...\n', flags);
	[status, out] = unix(sprintf('%s/tools/bowtie/bowtie %s %s %s > %s', ...
		ppath, flags, index_name, read_files{1}, alignments_file));
	if status ~= 0
		fprintf(1, '%s', out);
		error 'Bowtie read alignment failed.';
	end
	
	al_first = parse_alignments(alignments_file, out, bcols);
	al_first.ReadID = strcat(al_first.ReadID, '/1');
	
	fprintf(1, '-> Invoking Bowtie with flags "%s"...\n', flags);
	[status, out] = unix(sprintf('%s/tools/bowtie/bowtie %s %s %s > %s', ...
		ppath, flags, index_name, read_files{2}, alignments_file));
	if status ~= 0
		fprintf(1, '%s', out);
		error 'Bowtie read alignment failed.';
	end
	
	al_second = parse_alignments(alignments_file, out, bcols);
	al_second.ReadID = strcat(al_second.ReadID, '/2');
	
	al = struct;
	fields = fieldnames(al_first);
	for f = 1:length(fields)
		field = fields{f};
		eval(sprintf('al.%s = cat(1, al_first.%s, al_second.%s);', ...
			field, field, field));
	end
else
	% Single end read alignment.
	fprintf(1, '-> Invoking Bowtie with flags "%s"...\n', flags);
	
	[status, out] = unix(sprintf('%s/tools/bowtie/bowtie %s %s %s > %s', ...
		ppath, flags, index_name, read_files{1}, alignments_file));
	if status ~= 0
		fprintf(1, '%s', out);
		error 'Bowtie read alignment failed.';
	end
	
	if isempty(output_file)
		al = parse_alignments(alignments_file, out, bcols);
	end
end

alignments = al;







function al = parse_alignments(alignments_file, out, bcols)

al = struct;
al.TotalReads = NaN;
al.AlignedReads = NaN;

out = strread(out, '%s', 'delimiter', '\n');
for k = 1:length(out)
	match = regexp(out{k}, '^# reads processed: (\d+)', 'tokens');
	if length(match) == 1
		tmp = match{1}; al.TotalReads = str2double(tmp{1});
		continue;
	end
	
	match = regexp(out{k}, ...
		'# reads with at least one reported alignment: (\d+)', 'tokens');
	if length(match) == 1
		tmp = match{1}; al.AlignedReads = str2double(tmp{1});
		continue;
	end
end

if ~isnan(al.TotalReads) && ~isnan(al.AlignedReads)
	fprintf(1, '-> %d / %d (%.1f%%) reads with at least one alignment\n', ...
		al.AlignedReads, al.TotalReads, al.AlignedReads / al.TotalReads * 100);
end

parse_format = '';
file_col = 0;
if ismember(1, bcols), parse_format = [parse_format '%s']; end
if ismember(2, bcols), parse_format = [parse_format '%s']; end
if ismember(3, bcols), parse_format = [parse_format '%s']; end
if ismember(4, bcols), parse_format = [parse_format '%d']; end
if ismember(5, bcols), parse_format = [parse_format '%s']; end
if ismember(6, bcols), parse_format = [parse_format '%s']; end

fid = fopen(alignments_file);
data = textscan(fid, parse_format, 'Delimiter', '\t');
fclose(fid);

file_col = 0;
if ismember(1, bcols), file_col = file_col+1; al.ReadID = data{file_col}; end
if ismember(2, bcols)
	file_col = file_col + 1;
	strands = data{file_col};
	al.Strand = repmat(' ', length(strands), 1);
	for k = 1:length(strands), al.Strand(k) = strands{k}; end
end
if ismember(3, bcols), file_col = file_col+1; al.Target = data{file_col}; end
if ismember(4, bcols), file_col = file_col+1; al.Offset = data{file_col}; end
if ismember(5, bcols), file_col = file_col+1; al.Sequence = data{file_col}; end
if ismember(6, bcols), file_col = file_col+1; al.Quality = data{file_col}; end

