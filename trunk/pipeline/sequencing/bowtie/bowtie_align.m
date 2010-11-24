function [alignments_tmp, out] = bowtie_align(reads, index, options)

% If the index name is specified as a relative path, prefix it with the
% pipeline directory that holds all Bowtie indices.
if index(1) ~= '/'
	index = bowtie_index(index);
end

[flags, index_suffix] = bowtie_flags_for_reads(reads);

alignments_tmp = ptemp();
[status, out] = unix(sprintf('%s/tools/bowtie/bowtie %s %s %s%s %s > %s', ...
	ppath, options, flags, index, index_suffix, reads, alignments_tmp));

if status ~= 0
	fprintf(1, '%s', out);
	error 'Bowtie read alignment failed.';
end

