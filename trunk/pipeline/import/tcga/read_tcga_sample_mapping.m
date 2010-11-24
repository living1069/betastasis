function file_to_sample = read_tcga_sample_mapping(filepath)

mapping = fopen(filepath);
if mapping == -1
	fprintf(1, 'File %s could not be opened.', filepath);
	return;
end

fprintf(1, 'Reading TCGA sample-file mappings into memory...\n');
fprintf(1, 'Entries read: 0');

progress_len = 1;
entries = 0;

%sample_to_file = containers.Map();
file_to_sample = containers.Map();

while 1
	sample = '';
	
	while 1
		line = fgetl(mapping);
		if line == -1, break, end
		tokens = regexp(line, '<td.+>(TCGA-.+\d)</td>', 'tokens');
		if length(tokens) ~= 1, continue, end

		tokens = tokens{1}; sample = tokens{1};
		%fprintf(1, '%s\n', line);
		break;
	end
	
	tokens = [];
	no_files = 0;
	for k = 1:3
		line = fgetl(mapping);
		if line == -1, return, end

		tokens = regexp(line, '<td.+>(.+)</td>', 'tokens');
		if length(tokens) ~= 1
			%fprintf(1, 'Invalid file entry in TCGA metadata:\n', line);
			%fprintf(1, '%s\n', line);
			no_files = 1;
		end
	end
	
	if no_files, continue, end
	
	tokens = tokens{1};
	files = strrep(tokens{1}, '<br/>', '\n');
	files = strread(files, '%s', 'delimiter', '\n');
		
	entries = entries + 1;
	
	for k = 1:length(files)
		file = files{k};
		
		%if sample_to_file.isKey(sample)
		%	error 'Duplicate sample ID detected in TCGA metadata.';
		%end
		if file_to_sample.isKey(file)
			%fprintf(1, 'Duplicate file ID detected in TCGA metadata.\n');
			continue;
		end
		%sample_to_file(sample) = file;
		file_to_sample(file) = sample;
	end

	% Update the progress indicator that is shown to the user.
	if mod(entries, 100) == 0
		for j = 1:progress_len, fprintf(1, '\b'); end
		fprintf(1, '%d', entries);
		progress_len = length(int2str(entries));
	end
end

for j = 1:progress_len, fprintf(1, '\b'); end
fprintf(1, '%d\n', entries);

