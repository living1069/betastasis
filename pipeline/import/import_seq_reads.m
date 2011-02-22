
% IMPORT_SEQ_READS   Import high throughput sequence reads into the pipeline
%
%    IMPORT_SEQ_READS(DATASET) searches the current working directory 
%    high throughput sequencing read files, and imports any found read files
%    into the pipeline as a dataset. The name of the new dataset is
%    specified by the argument DATASET.
%
%    IMPORT_SEQ_READS(..., 'Recursive', TRUE) tells the function to also 
%    recursively look for read files in subfolders.

% Author: Matti Annala <matti.annala@tut.fi>

function [] = import_seq_reads(dataset, varargin)

tmp = ptemp;

recursive = false;
sequence_type = '';
platform = '';
paired_end = [];

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Platform')
		platform = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'SequenceType')
		sequence_type = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'PairedEnd')
		paired_end = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Recursive')
		recursive = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

reads = struct;
reads.Sample = struct;
reads.Sample.Filename = {};
reads.Sample.SequenceType = {};
reads.Resource = {};
reads.Type = 'Sequence reads';

if ~exist([ppath '/datasets'])
	error 'Cannot find pipeline datasets directory.';
end

single_files = {};
paired_files = {};

dirs = {'.'};

ds_dir = [ppath '/datasets/' flatten_str(dataset)];

while ~isempty(dirs)
	curr_dir = dirs{1};
	dirs = dirs(2:end);
	
	paired_prefixes = cell(0, 2);
	
	files = dir(curr_dir);
	
	for k = 1:length(files)
		if files(k).name(1) == '.', continue, end
			
		full_path = [curr_dir '/' files(k).name];
		
		if files(k).isdir
			if recursive, dirs{end+1} = full_path; end
			continue;
		end
			
		filename = files(k).name;
		
		if isempty(regexpi(filename, ...
			'\.(fa|fasta|csfasta|csfastq|fastq|fq|raw|sms|sam|bam)'))
			continue;   % Not a sequence file.
		end
		
		% Check if the file is paired end.
		tokens = regexpi(filename, '(.+)_(1|2)\.', 'tokens');
		if length(tokens) == 1 && (isempty(paired_end) || paired_end == true)
			token = tokens{1}; prefix = token{1}; side = str2double(token{2});
			
			match = find(strcmp(prefix, paired_prefixes(:, 1-(side-1)+1))); 
			if length(match) == 0
				paired_prefixes{end+1, side} = prefix;
				paired_files{end+1, side} = full_path;
			elseif length(match) == 1
				paired_prefixes{match, side} = prefix;
				paired_files{match, side} = full_path;
			end
		else
			single_files{end+1, 1} = full_path;
		end
	end
end
	




% Check that all paired files are accounted for.
for k = 1:length(paired_files)
	if isempty(paired_files{k, 1})
		fprintf(1, 'ERROR: Missing pair file for %s...\n', ...
			paired_files{k, 2});
		return;
	elseif isempty(paired_files{k, 2})
		fprintf(1, 'ERROR: Missing pair file for %s...\n', ...
			paired_files{k, 1});
		return;
	end
end





% Now we take all paired end files and add them to our dataset.
for k = 1:size(paired_files, 1)	
	
	first = path_strip_dir(paired_files{k, 1});
	second = path_strip_dir(paired_files{k, 2});
	
	if regexpi(first, '\.(fa|fasta|csfasta|csfastq|fastq|fq)')
		fprintf(1, 'Found paired end FASTA reads: %s and %s\n', first, second);
	elseif regexpi(filename, '\.raw')
		fprintf(1, 'Found raw reads: %s\n', filename);
		error 'Raw reads not supported yet.';
	elseif regexpi(filename, '\.sms')
		fprintf(1, 'Found Helicos SMS format reads: %s\n', filename);
		error 'SMS reads not supported yet.';
	elseif regexpi(filename, '\.sam')
		fprintf(1, 'Found SAM format alignments: %s\n', filename);
		error 'SAM reads not supported yet.';
	elseif regexpi(filename, '\.bam')
		fprintf(1, 'Found BAM format alignments: %s\n', filename);
		error 'BAM reads not supported yet.';
	else
		error 'Unrecognized file format.';
	end

	first_extracted = paired_files{k, 1};
	second_extracted = paired_files{k, 2};
	
	if regexpi(first, '\.gz$')
		fprintf(1, '-> Extracting %s...\n', first);
		first_extracted = [tmp '.1.extracted'];
		[status, ~] = unix(sprintf('gunzip -c %s > %s', ...
			paired_files{k, 1}, first_extracted));
		first = first(1:end-3);
	end
	
	if regexpi(second, '\.gz$')
		fprintf(1, '-> Extracting %s...\n', second);
		second_extracted = [tmp '.2.extracted'];
		[status, ~] = unix(sprintf('gunzip -c %s > %s', ...
			paired_files{k, 2}, second_extracted));
		second = second(1:end-3);
	end
	
	first_raw = [tmp '.1.raw'];
	second_raw = [tmp '.2.raw'];
	
	if regexpi(first, '\.(fa|fasta|csfasta|csfastq|fastq|fq)')
		
		fprintf(1, '-> Discarding headers and quality...\n');
		status = unix(sprintf([ppath ...
			'/sources/sequencing/transform/fasta_to_raw.py %s %s %s %s'], ...
			first_extracted, second_extracted, first_raw, second_raw));
		if status ~= 0
			error('fasta_to_raw.py returned an error:\n%s', out);
		end
	else
		error 'Unrecognized file format.';
	end
	
	if ~strcmp(first_extracted, paired_files{k, 1})
		safe_delete(first_extracted);
	end
	if ~strcmp(second_extracted, paired_files{k, 2})
		safe_delete(second_extracted);
	end
	
	[color, ~] = seq_read_type(first_raw);
	
	reads.Sample.Filename{end+1, 1} = first;
	reads.Resource{end+1, 1} = [flatten_str(dataset) '/' gen_uuid() '.gz'];

	if color
		reads.Sample.SequenceType{end+1, 1} = 'Colorspace, paired end';
	else
		reads.Sample.SequenceType{end+1, 1} = 'Nucleotide, paired end';
	end
	
	fprintf(1, '-> Storing reads in compressed format...\n');
	[~, ~] = mkdir(ds_dir);
	[status, ~] = unix(sprintf('gzip -c %s > %s', first_raw, ...
		[ppath '/datasets/' reads.Resource{end}(1:end-3) '.1' ...
		reads.Resource{end}(end-2:end)]));
	if status, error 'Gzip compression failed.'; end
	
	[status, ~] = unix(sprintf('gzip -c %s > %s', second_raw, ...
		[ppath '/datasets/' reads.Resource{end}(1:end-3) '.2' ...
		reads.Resource{end}(end-2:end)]));
	if status, error 'Gzip compression failed.'; end

	if ~strcmp(first_raw, first_extracted)
		safe_delete(first_raw);
	end
	if ~strcmp(second_raw, second_extracted)
		safe_delete(second_raw);
	end
end
	





	
% Now we take all single end files and add them to our dataset.
for k = 1:length(single_files)
	
	full_path = single_files{k};
	filename = path_strip_dir(full_path);
	
	extracted = full_path;
	
	if regexpi(filename, '\.gz$')
		extracted = [tmp '.extracted'];
		[status, ~] = unix(sprintf('gunzip -c %s > %s', full_path, extracted));
		filename = filename(1:end-3);
	end
	
	raw_file = [tmp '.raw'];
	
	if regexpi(filename, '\.(fa|fasta|csfasta|csfastq|fastq|fq)')
		fprintf(1, 'Found single end FASTA reads: %s\n', filename);
		
		fprintf(1, '-> Discarding headers and quality...\n');
		[status, out] = unix(sprintf([ppath ...
			'/sources/sequencing/transform/fasta_to_raw.py %s > %s'], ...
			extracted, raw_file));
		if status ~= 0
			error('fasta_to_raw.py returned an error:\n%s', out);
		end
		
	elseif regexpi(filename, '\.raw')
		fprintf(1, 'Found raw reads: %s\n', filename);
		error 'Raw reads not supported yet.';
		
	elseif regexpi(filename, '\.sam')
		fprintf(1, 'Found SAM format alignments: %s\n', filename);
		error 'SAM reads not supported yet.';
		
	elseif regexpi(filename, '\.sms')
		fprintf(1, 'Found Helicos SMS format reads: %s\n', filename);
		error 'SMS reads not supported yet.';
		
	elseif regexpi(filename, '\.bam')
		fprintf(1, 'Found BAM format alignments: %s\n', filename);
		error 'BAM reads not supported yet.';
		
	else
		error 'Unrecognized file format.';
	end
	
	if ~strcmp(extracted, full_path)
		safe_delete(extracted);
	end
	
	[color, ~] = seq_read_type(raw_file);
	
	reads.Sample.Filename{end+1, 1} = filename;
	reads.Resource{end+1, 1} = [flatten_str(dataset) '/' gen_uuid() '.gz'];

	if color
		reads.Sample.SequenceType{end+1, 1} = 'Colorspace, single end';
	else
		reads.Sample.SequenceType{end+1, 1} = 'Nucleotide, single end';
	end
	
	fprintf(1, '-> Storing reads in compressed format...\n');
	[~, ~] = mkdir(ds_dir);
	[status, ~] = unix(sprintf('gzip -c %s > %s', raw_file, ...
		[ppath '/datasets/' reads.Resource{end}]));
	if status, error 'Gzip compression failed.'; end
	
	if ~strcmp(raw_file, extracted)
		safe_delete(raw_file);
	end
end

S = size(paired_files, 1) + length(single_files);
if S == 0
	fprintf(1, 'No sequencing samples found.\n');
	return;
end

if ~isempty(platform)
	reads.Platform = repmat({platform}, S, 1);
end

metadata = reads;
save([ds_dir '/metadata.mat'], 'metadata');

