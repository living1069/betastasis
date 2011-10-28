function split = split_reads(reads, tag_len)

tmp = temporary('split_reads');

S = length(reads.url);

extracted = extract_reads(reads);
split = extracted;

read_files = seq_filenames(extracted);

for s = 1:S
	if ~rx(reads.paired{s}, 'single')
		error 'Reads are not single end.';
	end
	
	if ~rx(reads.format{s}, 'FASTA')
		error 'Only FASTA reads can be split at the moment.';
	end
	
	split.url{s} = [tmp reads.meta.sample_id{s}];
	
	for f = 1:length(read_files{s})
		
		left_split = [split.url{s} '_1.fa'];
		right_split = [split.url{s} '_2.fa'];
		
		[status, out] = unix(sprintf(['%s/sources/sequencing/transform/' ...
			'split_fasta_reads.py %s %d %s %s'], ...
			ppath, read_files{s}{f}, tag_len, left_split, right_split));
		if status ~= 0
			error('split_fasta_reads.py returned an error:\n%s', out);
		end
	end
	
	split.paired{s} = 'Paired';
end


