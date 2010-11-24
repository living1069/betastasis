function [flags, index_suffix] = bowtie_flags_for_reads(reads)

[color, quality] = seq_read_type(reads);

flags = '';
index_suffix = '';

if color == 1
	index_suffix = '_colorspace';
end

if color == 1 && quality == 1
	flags = '-C --solexa-quals --integer-quals';
elseif color == 1 && quality == 0
	flags = '-C -f';
elseif color == 0 && quality == 1
	flags = '';
elseif color == 0 && quality == 0
	flags = '-f';
end

