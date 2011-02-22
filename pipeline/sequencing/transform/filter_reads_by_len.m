function filtered = filter_reads_by_len(reads, len_range)

if length(len_range) == 1
	fprintf(1, '-> Discarding reads shorter than %d bases...\n', len_range(1));
	len_range(2) = 10000;
elseif length(len_range) == 2
	fprintf(1, ['-> Discarding reads with lengths outside range ' ...
			   '[%d, %d]...\n'], len_range(1), len_range(2));
else
	error 'Length range must be a vector of one or two elements.';
end

extracted = extract_reads(reads);

S = length(reads.Raw);

filtered = struct;
filtered.Meta = reads.Meta;
filtered.Raw = {};

for s = 1:S
	read_files = reads.Raw{s}.Paths;
	
	color = 0;
	if regexpi(reads.Meta.Sample.SequenceType, 'color')
		color = 1;
	end

	filtered.Raw{s} = FilePool;

	for f = 1:length(read_files)
	end
end
	


[color, quality] = seq_read_type(reads);

[status, output] = unix(sprintf( ...
	'%s/sources/sequencing/filter_reads_by_len.py %s %d %d > %s', ...
	ppath, reads, len_range(1) + 3 * color, len_range(2) + 3 * color, ...
	trimmed_tmp));
if status ~= 0
	error('Filtering of reads by length failed:\n%s\n', output);
end

