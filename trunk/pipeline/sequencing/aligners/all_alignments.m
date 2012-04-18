function al = all_alignments(alignments, varargin)

for al = iterate_alignments(alignments, 'Batch', Inf, varargin{:})
	return;
end

