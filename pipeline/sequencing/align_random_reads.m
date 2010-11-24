function al = align_random_reads(read_count, read_len, domain, index, varargin)

aligner = '';

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Aligner')
		aligner = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

N = read_count;
seq = randi([1 4], N, read_len);

if regexpi(domain, '^(nucleotide|tcga)$')
	char_map = [ 48 49 50 51 ];   % Colorspace ('0', '1', '2', '3')
	seq = [ (zeros(N, 1) + 84), char_map(seq) ];
elseif regexpi(domain, '^(color|colorspace|dinucleotide)$')
	char_map = [ 65 67 71 84 ];   % ACTG
	seq = char_map(seq);
else
	error 'Unsupported read domain.';
end

reads_tmp = ptemp;
write_seq_fasta(cellstr(char(seq)), reads_tmp);

fprintf(1, 'Aligning reads using Bowtie...\n');
al = align_reads(reads_tmp, index, varargin{:}, 'Columns', 'read,sequence');

