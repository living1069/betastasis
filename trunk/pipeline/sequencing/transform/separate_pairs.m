function paired = separate_pairs(reads, tag_len)

tmp = temporary('separate_pairs');

S = length(reads.url);

extracted = extract_reads(reads);
paired = extracted;

read_files = seq_filenames(extracted);

if any(~rx(reads.paired, 'single'))
	error 'All read samples must be single end.';
end

if any(~rx(reads.format, 'FASTA'))
	error 'Only FASTA reads can be split at the moment.';
end

paired.paired = repmat({'Paired'}, 1, S);

for s = 1:S
	paired.url{s} = [tmp reads.meta.sample_id{s}];
	
	for f = 1:length(read_files{s})
		left_split = [paired.url{s} '_1.fa'];
		right_split = [paired.url{s} '_2.fa'];
		
		[status, out] = unix(sprintf(['%s/sources/sequencing/transform/' ...
			'separate_fasta_pairs.py %s %s %s'], ...
			ppath, read_files{s}{f}, left_split, right_split));
		if status ~= 0
			error('separate_fasta_pairs.py returned an error:\n%s', out);
		end
	end
	
	split.paired{s} = 'Paired';
end


