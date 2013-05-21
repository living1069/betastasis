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
		data = read_uarray_sample_cel(file, probes);
	catch exception
		fprintf('Skipping sample due to import failure. Reason:\n');
		fprintf('%s\n', exception.message);
		continue;
	end

	if s == 1
		raw.mean = nan(length(data.mean), S);
	end
	
	raw.mean(:, s) = data.mean;
end

unix(sprintf('rm %s 2> /dev/null', decompressed));

