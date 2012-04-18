function trimmed = trim_reads(reads, trim_len)

use_pipes = true;
tmp = temporary('trim_reads');

S = length(reads.url);

fprintf('-> Trimming reads to a length of %d bases...\n', trim_len);

if any(~rx(reads.format, 'FASTA'))
	error 'Only FASTA reads can be trimmed at the moment.';
end

trimmed = reads;

for s = 1:S
	inputs = {};
	outputs = {};
	
	trimmed.url{s} = [tmp reads.meta.sample_id{s}];
	trimmed.format{s} = 'FASTA';
	
	if rx(reads.paired{s}, 'paired')
		inputs = strcat(reads.url{s}, {'_1.fa', '_2.fa'});
		outputs = strcat(trimmed.url{s}, {'_1.fa', '_2.fa'});
	else
		inputs = {[reads.url{s} '.fa']};
		outputs = {[trimmed.url{s} '.fa']};
	end
	
	if rx(reads.format{s}, 'gzip')
		for k = 1:length(inputs)
			inputs{k} = sprintf('<(gunzip -c %s)', [inputs{k} '.gz']);
		end
	end
	
	if use_pipes
		for k = 1:length(outputs)
			unix(['mkfifo ' outputs{k}]);
		end
	end
	
	color = rx(reads.space{s}, 'color');
	
	for k = 1:length(inputs)
		[status, out] = unix(sprintf( ...
			'%s/sources/sequencing/transform/trim_fasta_reads.py %s %d > %s &', ...
			ppath, inputs{k}, trim_len+1*color, outputs{k}));
		if status ~= 0
			error('trim_reads.py returned an error:\n%s', out);
		end
	end
end


