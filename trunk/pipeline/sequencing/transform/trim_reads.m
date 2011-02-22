function trimmed = trim_reads(reads, trim_len)

extracted = extract_reads(reads);

S = length(reads.Raw);

trimmed = struct;
trimmed.Meta = reads.Meta;
trimmed.Raw = {};

for s = 1:S
	read_files = reads.Raw{s}.Paths;
	
	color = 0;
	if regexpi(reads.Meta.Sample.SequenceType{s}, 'color')
		color = 1;
	end
	
	trimmed.Raw{s} = FilePool;
	
	for f = 1:length(read_files)
		trimmed_file = trimmed.Raw{s}.temp(sprintf('%d', f));
		[status, out] = unix(sprintf( ...
			'%s/sources/sequencing/transform/trim_raw_reads.py %s %d > %s', ...
			ppath, read_files{f}, trim_len+3*color, trimmed_file));
		if status ~= 0
			error('trim_reads.py returned an error:\n%s', out);
		end
	end
end


