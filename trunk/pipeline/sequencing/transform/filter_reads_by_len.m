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
	read_files = extracted.Raw{s}.Paths;
	
	color = ~isempty(regexpi(reads.Meta.Sequence.Space{s}, 'color'));

	filtered.Raw{s} = FilePool;

	for f = 1:length(read_files)
		filtered_file = filtered.Raw{s}.temp(sprintf('%d', f));
		[status, out] = unix(sprintf(['%s/sources/sequencing/transform/' ...
			'filter_raw_reads_by_len.py %s %d > %s'], ...
			ppath, read_files{f}, len_range(1)+1*color, ...
			len_range(2)+1*color, filtered_file));
		if status ~= 0, error('trim_reads.py returned an error:\n%s', out); end
	end
end

