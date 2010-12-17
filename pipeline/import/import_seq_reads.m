function [] = import_seq_reads(dataset, platform, varargin)

move_files = false;
recursive = false;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MoveFiles')
		move_files = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Recursive')
		recursive = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

reads = struct;
reads.Meta = struct;
reads.Meta.Type = 'Sequence reads';
reads.SequenceResource = {};

if ~exist([ppath '/datasets'])
	error 'Cannot find pipeline datasets directory.';
end

num_found = 0;
found_files = {};

dirs = {'.'};

while ~isempty(dirs)
	curr_dir = dirs{1};
	dirs = dirs(2:end);
	
	files = dir(curr_dir);
	for k = 1:length(files)
		if files(k).name(1) == '.', continue, end
			
		if files(k).isdir
			if recursive, dirs{end+1} = [curr_dir '/' files(k).name]; end
			continue;
		end
		
		tokens = regexpi(files(k).name, ...
			'(.+)\.(fa|fasta|csfasta|csfastq|fastq|fq|sms)', 'tokens');
		if isempty(tokens), continue, end
			
		filename = files(k).name;
		
		num_found = num_found + 1;
		found_files{num_found} = [curr_dir '/' filename];
		reads.Meta.Sample.Filename{num_found, 1} = filename;
		
		if regexpi(filename, '.+\.sms')
			fprintf(1, 'Found Helicos reads: %s\n', filename);
		else
			[color, ~] = seq_read_type([curr_dir '/' filename]);
				
			if color
				fprintf(1, 'Found colorspace reads: %s\n', filename);
			else
				fprintf(1, 'Found nucleotide space reads: %s\n', filename);
			end
		end
	end
end

if num_found == 0
	fprintf(1, 'No sequencing samples found.\n');
	return;
end

if length(unique(found_files)) ~= length(found_files)
	fprintf(1, 'WARNING: Found samples with an identical filename.\n');
end

reads.Meta.Platform = repmat({platform}, num_found, 1);

raw_dir = [ppath '/datasets/' flatten_str(dataset) '/raw'];
[~, ~] = mkdir(raw_dir);

for k = 1:length(found_files)
	if regexpi(found_files{k}, '.+\.sms')
		% For SMS files from the HeliScope platform, we first convert the
		% files into the standard FASTA format.
		setenv('HELICOS_TMPDIR', pipeline_config.TempDir);

		fprintf(1, 'Converting %s to FASTA...\n', ...
			path_strip_dir(found_files{k}));
		
		fasta_reads = [raw_dir '/' regexprep(path_strip_dir(found_files{k}), ...
			'^(.+).sms$', '$1', 'ignorecase')];
		[status, out] = unix(sprintf([ ...
			'%s/tools/helisphere/bin/sms2txt --fasta ' ...
			'--input_file %s --output_file_prefix %s'], ...
			ppath, found_files{k}, fasta_reads));
		if status ~= 0
			fprintf(1, '%s\n', out);
			error 'SMS to FASTA conversion failed.';
		end
		
		reads.SequenceResource{k} = ...
			[flatten_str(dataset) '/raw/' path_strip_dir(fasta_reads) '.fa'];
	else
		% For the standard FASTA format files, we either copy or move them into
		% the dataset repository.
		cmd = 'cp';
		if move_files, cmd = 'mv'; end
		[status, ~] = unix([cmd ' ' found_files{k} ' "' raw_dir '/"']);
		if status ~= 0
			error('Could not access file %s...', found_files{k});
		end
		
		reads.SequenceResource{k} = ...
			[flatten_str(dataset) '/raw/' found_files{k}];
	end
end

create_dataset(dataset, reads);

