function extracted = extract_reads(reads)

S = length(reads.Raw);

extracted = struct;
extracted.Meta = reads.Meta;
extracted.Raw = {};

for s = 1:S
	read_files = reads.Raw{s}.Paths;
	
	extracted.Raw{s} = FilePool;
	
	for f = 1:length(read_files)
		% Check if the reads even need to be extracted.
		if strcmpi(read_files{f}(end-2:end), '.gz')
			extracted_file = extracted.Raw{s}.temp(sprintf('%d', f));
			[status, out] = unix(sprintf('gunzip -c %s > %s', ...
				read_files{f}, extracted_file));
			if status ~= 0
				error('gunzip returned an error:\n%s', out);
			end
		else
			% FIXME: Maybe this should instead somehow inherit the
			% temp/static status.
			extracted.Raw{s}.static(read_files{f});
		end
	end
end


