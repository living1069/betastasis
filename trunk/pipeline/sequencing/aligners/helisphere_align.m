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

al = struct;
al.TotalReads = NaN;
al.FilteredReads = NaN;
al.AlignedReads = NaN;

flags = sprintf('--percent_error %d', 5 * max_mismatches);

if ~isempty(report_alignments) && ~isempty(allow_alignments) && ...
	report_alignments < allow_alignments
	fprintf(1, ['WARNING: ReportAlignments is higher than AllowAlignments.\n']);
	report_alignments = allow_alignments;
end

% Determine the type of index that the user wishes to align against.
if ischar(index)
	% If the index name is specified as a relative path, prefix it with the
	% pipeline directory that holds all Bowtie indices.
	if index(1) ~= '/'
		index = helisphere_index(index);
	else
		error 'Invalid Helisphere index.';
	end
else
	tmp_root = ptemp;
	fasta_tmp = [tmp_root '.fa'];
	write_seq_fasta(index, fasta_tmp);
	
	index = tmp_root;

	[status, ~] = unix(sprintf(['%s/tools/helisphere/bin/preprocessDB ' ...
		'--reference_file %s --out_prefix %s'], ppath, fasta_tmp, index));
	if status ~= 0, error 'Helisphere index construction failed.'; end
end

if isempty(output_file)
	alignments_file = tmp_pool.temp('align.bin');
else
	error 'Helisphere output capture not yet supported.';
end

filtered = tmp_pool.temp('.filtered.sms');

% FIXME: Make the reads.Raw{1}.Paths{1} less magic.
fprintf(1, '-> Discarding reads shorter than 24 bases...\n');
[status, out] = unix(sprintf(['%s/tools/helisphere/bin/filterSMS ' ...
	'--minlen 24 --input_file %s --output_file %s'], ...
	ppath, reads.Raw{1}.Paths{1}, filtered));
if status ~= 0, error('filterSMS failed:\n%s\n', out); end

% Capture the total number of reads from filterSMS output.
out = strread(out, '%s', 'delimiter', '\n');
for k = 1:length(out)
	tokens = regexp(out{k}, '(\d+) -> (\d+)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		
		if isnan(al.TotalReads), al.TotalReads = 0; end
		al.TotalReads = al.TotalReads + str2double(token{1});
		
		if isnan(al.FilteredReads), al.FilteredReads = 0; end
		al.FilteredReads = al.FilteredReads + str2double(token{2});
		
		fprintf(1, '-> %d / %d (%.1f%%) reads remain after filtering.\n', ...
			al.FilteredReads, al.TotalReads, ...
			al.FilteredReads / al.TotalReads * 100);
	end
end

fprintf(1, '-> Invoking Helisphere with flags "%s".\n', flags);
[status, out] = unix(sprintf(['%s/tools/helisphere/bin/indexDPgenomic %s ' ...
	'--read_file %s --reference_file %s.fa --data_base %s --output_file %s'],...
	ppath, flags, filtered, index, index, alignments_file));
if status ~= 0, error('indexDPgenomic failed:\n%s\n', out); end

filtered_alignments_file = strrep(alignments_file, '.bin', '.filtered.bin');
	
if ~isempty(allow_alignments)
	fprintf(1, '-> Discarding reads with too many alignments.\n');
	[status, out] = unix(sprintf(['%s/tools/helisphere/bin/filterAlign ' ...
		'--max_align %d --input_file %s --output_file %s'],...
		ppath, allow_alignments, alignments_file, filtered_alignments_file));
	if status ~= 0, error('filterAlign failed:\n%s\n', out); end
end

tab_al_tmp = regexprep(filtered_alignments_file, '\.bin$', '\.txt');

[status, out] = unix(sprintf(['%s/tools/helisphere/bin/align2txt ' ...
	'--full_reference --tab --input_file %s'], ppath, alignments_file));
if status ~= 0, error('align2txt failed:\n%s\n', out); end

fid = fopen(tab_al_tmp);
while 1
	line = fgetl(fid);
	if line == -1, error 'Alignment file ended abruptly.'; end
	if regexp(line, '^Reference_ID.*'), break, end
end

data = textscan(fid, '%s %s %d %d %*d %*d %*f %*d %*d %*d %*d %s %*s %s');
fclose(fid);

safe_delete(alignments_file);
safe_delete(filtered_alignments_file);
safe_delete(tab_al_tmp);

al.ReadID = data{2};
al.AlignedReads = length(unique(al.ReadID));  % FIXME: Slow?
fprintf(1, '-> Found alignments for %d (%.1f%%) reads.\n', ...
	al.AlignedReads, al.AlignedReads / al.TotalReads * 100);

strands = data{5};
al.Strand = repmat(' ', length(strands), 1);
for k = 1:length(strands), al.Strand(k) = strands{k}; end
	
al.Target = data{1};
al.Offset = data{3};

% FIXME: We should probably differentiate between the read sequence and
% reference sequence.
al.Sequence = strrep(data{6}, '-', '');

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

