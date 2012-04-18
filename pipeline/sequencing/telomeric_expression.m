function telomere_expr = telomeric_expression(reads, varargin)

anchor_len = NaN;
%anchor_len = 20;
max_mismatches = 2;

S = length(reads.url);

human_telomere = repmat('TTAGGG', 1, 20);

telomere_expr = nan(1, S);

for s = 1:S
	fprintf(1, 'Analyzing sample %s for telomeric expression...\n', ...
		reads.meta.sample_id{s});
	
	if anchor_len > 0
		fprintf(1, 'Splitting reads into %d bp start and end tags...\n', ...
			anchor_len);
		split = split_reads(filter(reads, s), anchor_len);
	else
		split = filter(reads, s);
	end

	alignments = bowtie_align(split, {human_telomere}, ...
		sprintf('-v%d', max_mismatches));
		
	al = all_alignments(alignments);
	for k = 1:length(al.read), al.read{k} = al.read{k}(1:end-2); end
	read_num = str2double(al.read);
	
	telomeric_reads = false(1, max(read_num));
	telomeric_reads(read_num) = true;
	
	telomere_expr(s) = sum(telomeric_reads);
end


