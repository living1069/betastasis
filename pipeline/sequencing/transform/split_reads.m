function split = split_reads(reads, tag_len)

extracted = extract_reads(reads);

S = length(reads.Raw);

split = struct;
split.Meta = reads.Meta;
split.Raw = {};

for s = 1:S
	read_files = extracted.Raw{s}.Paths;
	
	if isempty(regexpi(extracted.Meta.Sequence.Paired, 'single')) || ...
		length(read_files) ~= 1
		error 'Reads do not look single end.';
	end
	
	split.Raw{s} = FilePool;
	
	for f = 1:length(read_files)
		left_split = split.Raw{s}.temp('1');
		right_split = split.Raw{s}.temp('2');
		
		[status, out] = unix(sprintf( ...
			'%s/sources/sequencing/transform/split_raw_reads.py %s %d %s %s',...
			ppath, read_files{f}, tag_len, left_split, right_split));
		if status ~= 0
			error('split_raw_reads.py returned an error:\n%s', out);
		end
	end
	
	split.Meta.Sequence.Paired{s} = 'Paired end';
end


