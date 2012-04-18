function paths = access_reads(reads, varargin)

if length(reads.url) ~= 1
	error 'access_reads() should be called with one read at a time.';
end

s = 1;

if rx(reads.paired{s}, 'paired')
	paths{1} = [reads.url{s} '_1'];
	paths{2} = [reads.url{s} '_2'];
else
	paths{1} = reads.url{s};
end

if rx(reads.format{s}, 'FASTA')
	for k = 1:length(paths), paths{k} = [paths{k} '.fa']; end
elseif rx(reads.format{s}, 'FASTQ')
	for k = 1:length(paths), paths{k} = [paths{k} '.fq']; end
end

if rx(reads.format{s}, 'gzip')
	for k = 1:length(paths), paths{k} = ['<(gunzip -c ' paths{k} '.gz)']; end
end

