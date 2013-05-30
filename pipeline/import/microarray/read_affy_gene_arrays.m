function raw = read_affy_gene_arrays(sample_files, probes_path)

if ~iscellstr(sample_files)
	sample_files = find_files(sample_files);
end

S = length(sample_files);

raw = struct;
raw.rows = struct;
raw.meta.sample_id = sample_files(:);

fprintf('Loading probes...\n');
probes = load(probes_path); probes = probes.probes;

decompressed = [temporary('read_affy_gene_arrays') 'decompressed.CEL'];

for s = 1:S
	fprintf('Importing %s...\n', sample_files{s});
	
	file = sample_files{s};
	if rx(file, '^(.+)\.gz$')
		unix(sprintf('gunzip -c %s > %s', file, decompressed));
		file = decompressed;
	end
	
	try
		cel = affyread(file);
	catch exception
		fprintf('Skipping sample due to import failure. Reason:\n%s\n', ...
			exception.message);
		continue;
	end
	
	data = zeros(length(probes.sequence), 1, 'single');
	idx = probes.ypos * cel.Cols + probes.xpos + 1;
	data = cel.Probes(idx, 3);

	if s == 1
		raw.mean = nan(length(data), S);
		raw.rows.gene = probes.target;
	end
	raw.mean(:, s) = data;
end

unix(sprintf('rm %s 2> /dev/null', decompressed));

