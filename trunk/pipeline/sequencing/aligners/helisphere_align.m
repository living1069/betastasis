function al = helisphere_align(reads, index, varargin)

global pipeline_config;

tmp_pool = FilePool;

max_mismatches = 0;
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

flags = sprintf('--percent_error %d', 5 * max_mismatches);

if ~isempty(report_alignments) && ~isempty(allow_alignments)
	fprintf(1, ['WARNING: Specified AllowAlignments and ReportAlignments. ' ...
		'Using only the AllowAlignments option.\n']);
	report_alignments = [];
end

% If the index name is specified as a relative path, prefix it with the
% pipeline directory that holds all Helisphere indices.
if index(1) ~= '/'
	index = helisphere_index(index);
end

if isempty(output_file)
	alignments_file = tmp_pool.temp('align.bin');
else
	error 'Helisphere output capture not yet supported.';
end

filtered = tmp_pool.temp('.filtered.sms');

% FIXME: Make the reads.Raw{1}.Paths{1} less magic.
fprintf(1, '-> Discarding reads shorter than 20 bases...\n');
[status, out] = unix(sprintf(['%s/tools/helisphere/bin/filterSMS ' ...
	'--minlen 20 --input_file %s --output_file %s'], ...
	ppath, reads.Raw{1}.Paths{1}, filtered));
if status ~= 0, error('filterSMS failed:\n%s\n', out); end

fprintf(1, '-> Invoking Helisphere with flags "%s".\n', flags);
status = unix(sprintf(['%s/tools/helisphere/bin/indexDPgenomic %s ' ...
	'--read_file %s --reference_file %s.fa --data_base %s --output_file %s'],...
	ppath, flags, filtered, index, index, alignments_file));
if status ~= 0, error('Helisphere read alignment failed:\n%s\n', out); end

tab_al_tmp = regexprep(alignments_file, '\.bin$', '\.txt')

[status, out] = unix(sprintf(['%s/tools/helisphere/bin/align2txt ' ...
	'--full_reference --tab --input_file %s'], ppath, alignments_file));
if status ~= 0, error('align2txt failed:\n%s\n', out); end

fid = fopen(tab_al_tmp);
while 1
	line = fgetl(fid);
	if line == -1, error 'Alignment file ended abruptly.'; end
	if regexp(line, '^Reference_ID.*'), break, end
end

data = textscan(fid, '%s %s %d %d %*d %*d %*f %*d %*d %*d %*d %s %*s %*s');
fclose(fid);

al = struct;
if ismember(1, bcols)
	al.ReadID = data{2};
elseif ismember(2, bcols)
	strands = data{5};
	al.Strand = repmat(' ', length(strands), 1);
	for k = 1:length(strands), al.Strand(k) = strands{k}; end
elseif ismember(3, bcols)
	al.Target = data{1};
elseif ismember(4, bcols)
	al.Offset = data{3};
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

index_path = [ppath '/tools/helisphere/indexes/' flatten_str(organism.Name) ...
	'/' flatten_str(organism.Version) '/' index_path];

