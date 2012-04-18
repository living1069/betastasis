function iter = iterate_alignments(al, varargin)

if length(al.url) > 1
	error 'iterate_alignments() only works on singular alignment samples.';
end

if rx(al.paired, 'Paired'), error 'Only single end reads supported.'; end

s = 1;

if rx(al.format{s}, 'Bowtie')
	iter = BowtieAlignments([al.url{s} '.bowtie'], al.total_alignments(s), ...
		varargin{:});
elseif rx(al.format{s}, 'BAM')
	error 'BAM alignments not supported yet.';
else
	error 'Unsupported alignment format.';
end

