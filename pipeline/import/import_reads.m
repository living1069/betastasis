
% IMPORT_READS   Import high throughput sequence reads into the pipeline
%
%    IMPORT_READS(DATASET) searches the current working directory 
%    high throughput sequencing read files, and imports any found read files
%    into the pipeline as a dataset. The name of the new dataset is
%    specified by the argument DATASET.
%
%    IMPORT_READS(..., 'Recursive', TRUE) tells the function to also 
%    recursively look for read files in subfolders.

% Author: Matti Annala <matti.annala@tut.fi>

function reads = import_reads(varargin)

tmp = temporary('reads');

search_dir = '.';
recursive = false;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'recurs')
		recursive = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

reads.meta.type = 'Sequence reads';

single_files = {};
paired_prefixes = cell(0, 2);
paired_files = cell(0, 2);

dirs = {search_dir};

while ~isempty(dirs)
	curr_dir = dirs{1};
	dirs = dirs(2:end);
	
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
			'\.(fa|fasta|csfasta|csfastq|fastq|fq|sms|sam|bam)(\.gz)?$'))
			continue;   % Not a sequence file.
		end
		
		% Check if the file is paired end.
		tokens = regexpi(filename, '(.+)_(1|2)\.', 'tokens');
		if length(tokens) == 1
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
keep = true(size(paired_files, 1), 1);
for k = 1:size(paired_files, 1)
	if isempty(paired_files{k, 1})
		fprintf(1, 'WARNING: Missing pair file for %s...\n', ...
			paired_files{k, 2});
		single_files{end+1, 1} = paired_files{k, 2};
		keep(k) = false;
	elseif isempty(paired_files{k, 2})
		fprintf(1, 'WARNING: Missing pair file for %s...\n', ...
			paired_files{k, 1});
		single_files{end+1, 1} = paired_files{k, 1};
		keep(k) = false;
	end
end
paired_files = paired_files(keep, :);




S = 0;   % Number of sequence samples found.

% Now we take all paired end files and add them to our dataset.
for k = 1:size(paired_files, 1)	
	
	first_full = absolutepath(paired_files{k, 1});
	first = path_strip_dir(paired_files{k, 1});
	
	paired_end_rx = '^(.+)_.\.(fa|fasta|csfasta|fq|fastq)($|\..+)';
	if ~rx(first, paired_end_rx)
		fprintf(1, 'WARNING: Unrecognized file format "%s".\n', first);
		continue;
	end
	
	S = S + 1;
	[reads.format{S}, reads.space{S}] = decipher_file_format(first);
	reads.meta.sample_id{S} = regexprep(first, paired_end_rx, '$1');
	reads.paired{S} = 'Paired';
	reads.url{S} = regexprep(first_full, paired_end_rx, '$1');
	
	if rx(first, '\.(gz|gzip)$')
		reads.format{S} = [reads.format{S} ' (gzip)'];
	end
	
	fprintf(1, 'Found paired end %s reads: %s\n', reads.format{S}, ...
		regexprep(first, paired_end_rx, '$1_*.$2$3'));
end







% Now we take all single end files and add them to our dataset.
for k = 1:length(single_files)
	
	full_path = absolutepath(single_files{k});
	filename = path_strip_dir(full_path);
	
	single_end_rx = '^(.+)\.(fa|fasta|csfasta|fq|fastq)($|\.).*';
	if ~rx(filename, single_end_rx)
		fprintf(1, 'WARNING: Unrecognized file format "%s".\n', filename);
		continue;
	end
	
	S = S + 1;
	[reads.format{S}, reads.space{S}] = decipher_file_format(filename);
	reads.meta.sample_id{S} = regexprep(filename, single_end_rx, '$1');
	reads.paired{S} = 'Single';
	reads.url{S} = regexprep(full_path, single_end_rx, '$1');
	
	if rx(filename, '\.(gz|gzip)$')
		reads.format{S} = [reads.format{S} ' (gzip)'];
	end
	
	fprintf(1, 'Found single end %s reads: %s\n', reads.format{S}, filename);
end



if S == 0
	fprintf(1, 'No sequencing samples found.\n');
	return;
end






function [format, space] = decipher_file_format(filename)

if regexpi(filename, '\.(fa|fasta)')
	format = 'FASTA'; space = 'Nucleotide';
elseif regexpi(filename, '\.csfasta')
	format = 'FASTA'; space = 'Color';
elseif regexpi(filename, '\.(fq|fastq)')
	format = 'FASTQ'; space = 'Nucleotide';
elseif regexpi(filename, '\.csfastq')
	format = 'FASTQ'; space = 'Color';
else
	error 'Invalid file format after regex.';
end


