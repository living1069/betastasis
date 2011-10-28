function extracted = extract_reads(reads)

S = length(reads.url);

tmp = temporary('extract_reads');
extracted = reads;

read_files = seq_filenames(reads);

for s = 1:S
	if rx(reads.format{s}, 'gzip')
		extracted.format{s} = strrep(extracted.format{s}, ' (gzip)', '');
		extracted.url{s} = [tmp reads.meta.sample_id{s}];
		pipe_files = seq_filenames(filter(extracted, s));
		for f = 1:length(read_files{s})
			decompress_pipe(read_files{s}{f}, pipe_files{1}{f});
		end
	end
end


