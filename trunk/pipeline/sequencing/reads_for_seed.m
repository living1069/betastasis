function [] = reads_for_seed(reads, seed)

threshold = 10;

seed = upper(seed);
if isempty(regexpi(seed, '[ACGT]+'))
	error 'Seed must be given as a string of nucleotides.';
end

% Write all K-mer sequences to a FASTA file.
kmers = dec2base(0:4^6-1, 4);
kmers(kmers == '0') = 'A';
kmers(kmers == '1') = 'C';
kmers(kmers == '2') = 'G';
kmers(kmers == '3') = 'T';

seed_ext = [repmat(seed, 4^6, 1) kmers];



for s = 1:length(reads.Raw)
	trimmed = trim_reads(filter_query(reads, s), length(seed) + 6);
	
	al = align_reads(trimmed, seed_ext, 'MaxMismatches', 0, ...
		'Columns', 'target');
	
	num_reads = zeros(size(kmers, 1), 1),
	targets = str2double(al.Target);
	
	for t = 1:length(targets)
		num_reads(t) = num_reads(t) + 1;
	end
	
	fprintf(1, 'Most prevalent reads for sample %d:\n', s);
	
	good = find(num_reads >= threshold);
	[~, order] = sort(num_reads(good), 'descend');
	
	for k = 1:length(order)
		idx = good(order(k));
		fprintf(1, '%s|%s: %d reads\n', seed, kmers(idx, :), num_reads(idx));
	end
	
	fprintf(1, '\n');
end




