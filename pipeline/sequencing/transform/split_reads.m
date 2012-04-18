function split = split_reads(reads, anchor_len)

use_pipes = true;

tmp = temporary('split_reads');

S = length(reads.url);

split = reads;

if any(~rx(reads.format, 'FASTA'))
	error 'Only FASTA reads can be split at the moment.';
end

for s = 1:S
	split.url{s} = [tmp reads.meta.sample_id{s}];
	split.paired{s} = 'Paired';
	split.format{s} = 'FASTA';
	
	left_split = [split.url{s} '_1.fa'];
	right_split = [split.url{s} '_2.fa'];
	
	if use_pipes
		unix(['mkfifo ' left_split]);
		unix(['mkfifo ' right_split]);
	end
	
	cmd = 'python ~/pipeline/sequencing/transform/trim_fasta_reads.py';
	
	if rx(reads.format{s}, 'gzip')
		input_format = '<(gunzip -c %s)';
	else
		input_format = '%s';
	end
	
	if rx(reads.paired{s}, 'paired')
		input_1 = sprintf(input_format, [reads.url{s} '_1.fa']);
		input_2 = sprintf(input_format, [reads.url{s} '_2.fa']);
		
		[status, out] = unix(sprintf( ...
			'cat <(%s %s %d /1) <(%s %s %d /1) > %s &', ...
			cmd, input_1, anchor_len, cmd, input_2, anchor_len, left_split));
		if status ~= 0, error('trim_reads.py returned an error:\n%s', out); end
			
		[status, out] = unix(sprintf( ...
			'cat <(%s %s %d /2) <(%s %s %d /2) > %s &', ...
			cmd, input_1, -anchor_len, cmd, input_2, -anchor_len, right_split));
		if status ~= 0, error('trim_reads.py returned an error:\n%s', out); end
			
	else
		input = sprintf(input_format, [reads.url{s} '.fa']);
		
		[status, out] = unix(sprintf('%s %s %d /1 > %s &', ...
			cmd, input, anchor_len, left_split));
		if status ~= 0, error('trim_reads.py returned an error:\n%s', out); end
		
		[status, out] = unix(sprintf('%s %s %d /2 > %s &', ...
			cmd, input, -anchor_len, right_split));
		if status ~= 0, error('trim_reads.py returned an error:\n%s', out); end
	end
end


