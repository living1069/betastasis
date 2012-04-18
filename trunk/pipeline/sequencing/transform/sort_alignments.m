function sorted = sort_alignments(alignments, varargin)

in_place = true;

tmp = temporary('sorted_alignments');

if any(~rx(alignments.format, 'BAM'))
	error 'Only BAM format alignments can currently be sorted.';
end

sorted = alignments;

for s = 1:length(alignments.url)
	if ~in_place
		sorted.url{s} = [tmp alignments.meta.sample_id{s}];
	end
	
	[status, out] = unix(sprintf('samtools sort %s.bam %s', ...
		alignments.url{s}, sorted.url{s}));
	if status ~= 0, error('samtools sort returned an error:\n%s', out); end
		
	[status, out] = unix(sprintf('samtools index %s.bam', alignments.url{s}));
	if status ~= 0, error('samtools index returned an error:\n%s', out); end
end



