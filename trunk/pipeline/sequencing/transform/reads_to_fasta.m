function raw = reads_to_raw(reads)

extracted = extract_reads(reads);
S = length(reads.Raw);

raw_file = ptemp;

raw = struct;
raw.Meta = reads.Meta;
raw.Raw = {};

for s = 1:S
	fprintf(1, 'Converting sample %s to raw format:\n', ...
		reads.Meta.Sample.Filename{s});
	
	read_files = extracted.Raw{s}.Paths;
	raw.Raw{s} = FilePool;
	
	for f = 1:length(read_files)
		compressed = raw.Raw{s}.temp(sprintf('%d.gz', f));
	
		if strcmpi(reads.Meta.Sequence.Format{s}, 'FASTA')
			fprintf(1, '-> Discarding FASTA headers and quality...\n');
			[status, out] = unix(sprintf([ppath ...
				'/sources/sequencing/transform/fasta_to_raw.py %s > %s'], ...
				read_files{f}, raw_file));
			if status, error('fasta_to_raw.py returned an error:\n%s', out); end
		else
			error('Cannot convert reads of type %s to raw.', ...
				reads.Meta.Sequence.Format{s});
		end
		
		fprintf(1, '-> Storing reads in compressed format...\n');
		[status, out] = unix(sprintf('gzip -c %s > %s', raw_file, compressed));
		if status, error('Gzip compression failed:\n%s\n.', out); end
	end
	
	reads.Meta.Sequence.Format{s, 1} = 'Raw';
end


