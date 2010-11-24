function results = helisphere_align(reads, index, varargin)

global pipeline_config;

max_mismatches = 0;
read_len_filter = [];
trim_len = NaN;
report_alignments = [];
allow_alignments = [];
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
		fprintf(1, ['WARNING: ReportAlignments option not supported ' ...
			'by Helisphere aligner.\n']);
		report_alignments = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'AllowAlignments')
		fprintf(1, ['WARNING: AllowAlignments option not supported ' ...
			'by Helisphere aligner.\n']);
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

setenv('HELICOS_TMPDIR', pipeline_config.TempDir);
pwd = cd;
cd(pipeline_config.TempDir);

if regexpi(reads, '.*\.sms')
	fasta_reads = sms_to_fasta(reads);
else
	fasta_reads = reads;
end

[color, quality] = seq_read_type(fasta_reads);
if color, error 'Colorspace reads are not supported by Helisphere.'; end
if quality, error 'Reads with quality are not supported by Helisphere.'; end

flags = '--input_file_type fasta';
flags = sprintf('%s --percent_error %d', flags, 5 * max_mismatches);

if ~isempty(report_alignments) && ~isempty(allow_alignments)
	fprintf(1, ['WARNING: Specified AllowAlignments and ReportAlignments. ' ...
		'Using only the AllowAlignments option.\n']);
	report_alignments = [];
end

if ~isempty(read_len_filter)
	filtered_tmp = filter_reads_by_len(fasta_reads, read_len_filter);
else
	filtered_tmp = fasta_reads;
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

% If the index name is specified as a relative path, prefix it with the
% pipeline directory that holds all Helisphere indices.
if index(1) ~= '/'
	index = helisphere_index(index);
end

fprintf(1, 'Using flags "%s" when invoking Helisphere.\n', flags);

if isempty(output_file)
	alignments_file = [tmp '.bin'];
else
	if ~strcmp(output_file(end-3:end), '.txt')
		fprintf(1, ['WARNING: Helisphere output file must have .txt ' ...
		            'extension. Appending it now...']);
		alignments_file = [alignments_file '.txt'];
	end

	alignments_file = regexprep(output_file, '\.txt$', '.bin');
end

[status, out] = unix(sprintf(['%s/tools/helisphere/bin/indexDPgenomic %s ' ...
	'--read_file %s --reference_file %s.fa --data_base %s --output_file %s'],...
	ppath, flags, trimmed_reads_tmp, index, index, alignments_file));
if status ~= 0
	fprintf(1, '%s', out);
	error 'Helisphere read alignment failed.';
end

if ~isempty(read_len_filter)
	safe_delete(filtered_tmp);
end
	
if trim_len > 0
	safe_delete(trimmed_reads_tmp);
end

tab_al_tmp = regexprep(alignments_file, '\.bin$', '\.txt')

status = unix(sprintf(['%s/tools/helisphere/bin/align2txt --tab ' ...
	'--input_file %s'], ppath, alignments_file));
if status ~= 0
	error 'Could not convert binary alignments to ASCII.';
end

safe_delete(alignments_file);

% If the user wants to keep the output file, then we assume that he is not
% interested in the parsed file contents.
if ~isempty(output_file)
	fprintf(1, 'WARNING: Raw Helisphere output file requested.\n');
	return;
end

fid = fopen(tab_al_tmp);
while 1
	line = fgetl(fid);
	if line == -1, error 'Alignment file ended abruptly.'; end
	if regexp(line, '^Reference_ID.*'), continue, end
end

data = textscan(fid, '%d %d %d %d %*d %*d %*f %*d %*d %*d %*d %s %*s %*s');
fclose(fid);

safe_delete(tab_al_tmp);

results = struct;
if ismember(1, bcols)
	results.ReadID = cell(size(data{2}));
	for k = 1:length(results.ReadID)
		results.ReadID{k} = num2str(data{2}(k));
	end
elseif ismember(2, bcols)
	results.Strand = data{5};
elseif ismember(3, bcols)
	results.Target = organism.Transcripts.Name(data{1});
elseif ismember(4, bcols)
	results.Offset = data{3};
elseif ismember(5, bcols)
	error 'Sequence information requested from Helisphere.';
elseif ismember(6, bcols)
	error 'Quality information requested from Helisphere.';
elseif ismember(8, bcols)
	error 'List of alignment mismatches requested from Helisphere.';
end

cd(pwd);














function index_path = helisphere_index(feature)

global organism;

feature = lower(feature);
if regexp(feature, '^(genome|transcripts|exons|mirnas|pre_mirnas)$')
	index_path = feature;
else
	error('Helisphere index requested for unsupported feature "%s".', feature);
end

index_path = [ppath '/tools/helisphere/indexes/' ...
	flatten_str(organism.Version) '/' index_path];

