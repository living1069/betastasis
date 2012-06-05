function sorted = sort_alignments(alignments, varargin)

in_place = true;

tmp = temporary('sorted_alignments');

if any(~rx(alignments.format, 'BAM'))
	error 'Only BAM format alignments can be sorted.';
end

sorted = alignments;

for s = 1:length(alignments.url)
	if ~in_place
		sorted.url{s} = [tmp alignments.meta.sample_id{s}];
	end
	
	filename = basename(alignments.url{s});
	
	% Check if the BAM file is already sorted.
	[status, out] = unix(sprintf('samtools view -H %s.bam', alignments.url{s}));
	if status == 0 && rx(out, 'SO:coordinate')
		fprintf('%s.bam is already sorted.\n', filename);
	else
		fprintf('Sorting %s...\n', filename);
		[status, out] = unix(sprintf('samtools sort %s.bam %s', ...
			alignments.url{s}, sorted.url{s}));
		if status ~= 0, error('samtools sort returned an error:\n%s', out); end
	end
	
	% Check if the BAM file is already indexed.
	if exist([alignments.url{s} '.bam.bai']) || ...
		exist([alignments.url{s} '.bai'])
		fprintf('%s.bam already has an index.\n', filename);
	else
		fprintf('Building index for %s...\n', filename);
		[status, out] = unix(sprintf('samtools index %s.bam', ...
			alignments.url{s}));
		if status ~= 0, error('samtools index returned an error:\n%s', out); end
	end
end



