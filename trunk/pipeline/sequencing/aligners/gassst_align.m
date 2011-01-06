function results = gassst_align(reads, index, varargin)

global pipeline_config;

max_mismatches = 0;
read_len_filter = [];
trim_len = NaN;
report_alignments = [];
allow_alignments = [];
max_threads = pipeline_config.MaxThreads;
columns = 'read,strand,ref,offset';
output_file = '';

tmp = ptemp;

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
		continue;
	end
	
	if strcmpi(varargin{k}, 'ReportAlignments')
		report_alignments = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'AllowAlignments')
		fprintf(1, ['WARNING: AllowAlignments option not supported ' ...
			'by GASSST aligner.\n']);
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

[color, quality] = seq_read_type(reads);
if color, error 'Colorspace reads are not supported by GASSST.'; end
if quality, error 'Reads with quality are not supported by GASSST.'; end
	
if ischar(index)
	% If the index name is specified as a relative path, prefix it with the
	% pipeline directory that holds all GASSST indices.
	if index(1) ~= '/'
		index = gassst_index(index);
	end
	index_name = index;
else
	index_name = [ptemp '.fa'];
	write_seq_fasta(index, index_name);
end

flags = '';
flags = sprintf('%s -p %d', flags, 100 - 5 * max_mismatches);
flags = sprintf('%s -n %d', flags, max_threads);

if ~isempty(report_alignments) && ~isempty(allow_alignments)
	fprintf(1, ['WARNING: Specified AllowAlignments and ReportAlignments. ' ...
		'Using only the AllowAlignments option.\n']);
	report_alignments = [];
end

if ~isempty(report_alignments)
	flags = sprintf('%s -h %d', flags, report_alignments);
end

if ~isempty(read_len_filter)
	filtered_tmp = filter_reads_by_len(reads, read_len_filter);
else
	filtered_tmp = reads;
end

if trim_len > 0
	trimmed_reads_tmp = [tmp '.trimmed'];
	
	[color, quality] = seq_read_type(reads);
	
	fprintf(1, 'Trimming reads to a length of %d bp...\n', trim_len);
	[status, ~] = unix(sprintf( ...
		'%s/sources/sequencing/trim_reads.py %s %d > %s', ...
		ppath, filtered_tmp, trim_len + 3 * color, trimmed_reads_tmp));
else
	trimmed_reads_tmp = filtered_tmp;
end

fprintf(1, 'Invoking GASST with flags "%s".\n', flags);

if isempty(output_file)
	alignments_file = [tmp '.aligned'];
else
	alignments_file = output_file;
end

[status, out] = unix(sprintf( ...
	'%s/tools/gassst/Gassst %s -i %s -d %s -o %s', ...
	ppath, flags, trimmed_reads_tmp, index_name, alignments_file));
if status ~= 0
	fprintf(1, '%s', out);
	error 'GASSST read alignment failed.';
end

results = struct;

out = strread(out, '%s', 'delimiter', '\n');
for k = 1:length(out)
	tokens = regexp(out{k}, [regexptranslate('escape', trimmed_reads_tmp) ...
		':\s*(\d+) sequences'], 'tokens');
	if length(tokens) ~= 1, continue, end
	
	token = tokens{1};
	results.TotalReads = str2double(token{1});
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
	fprintf(1, 'WARNING: Raw GASSST output file requested.\n');
	return;
end

fid = fopen(alignments_file);
while 1
	line = fgetl(fid);
	if length(line) < 3
		error 'Alignment file ended abruptly.';
	end
	if ~strcmp(line(1:3), '@SQ'), break, end
end

if fseek(fid, -(length(line)+1), 'cof') ~= 0
	error 'fseek() to the alignment section failed.';
end

data = textscan(fid, '%s %s %d %d %s %*d %d %s %*d');
fclose(fid);

results.ReadID = data{1};
results.Strand = repmat('+', length(data{3}), 1);
results.Strand(data{3} == 1) = '-';
results.Target = data{2};
results.Offset = data{4};
results.Sequence = data{7};
results.ErrorCount = data{6};
results.Mismatches = data{5};

safe_delete(alignments_file);











function index_path = gassst_index(feature)

global organism;

feature = lower(feature);
if regexp(feature, '^(genome|transcripts|exons|mirnas|pre_mirnas)$')
	index_path = feature;
else
	error('GASSST index requested for unsupported feature "%s".', feature);
end

index_path = [ppath '/tools/gassst/indices/' flatten_str(organism.Name) '/' ...
	flatten_str(organism.Version) '/' index_path '.fa'];

