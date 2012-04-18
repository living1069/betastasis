
function microbes = microbial_expression(reads)

max_genome_mismatches = 1;

S = length(reads.url);
M = 0;

% We assume that transcriptome aligned reads have already been filtered out,
% or that we are dealing with DNA-seq and they need not be filtered out.
fprintf(1, 'Discarding reads that align to the genome...\n');
[~, unaligned] = bowtie_align(reads, 'genome', ...
	sprintf('-v%d', max_genome_mismatches), 'Output', false);

alignments = bowtie_align(unaligned, ...
	'/home/annalam/organisms/microbes/microbes', '-v0 -k100');
	
microbe_to_idx = containers.Map;

microbes.meta = reads.meta;
microbes.expr = zeros(0, S);
microbes.read_seqs = cell(0, S);

for s = 1:S
	fprintf(1, 'Calculating microbial expression in sample %s...\n', ...
		reads.meta.sample_id{s});
	
	for al = iterate_alignments(filter(alignments, s))
		for k = 1:length(al.target)
			m = al.target{k};
			if ~microbe_to_idx.isKey(m)
				M = M + 1; microbe_to_idx(m) = M;
				microbes.name{M, 1} = m;
				microbes.expr(M, s) = 0;
				microbes.read_seqs{M, s} = [];
			end
			m = microbe_to_idx(m);
			
			microbes.expr(m, s) = microbes.expr(m, s) + 1;
			if isempty(microbes.read_seqs{m, s})
				microbes.read_seqs{m, s} = {};
			end
			microbes.read_seqs{m, s}{end+1, 1} = al.sequence{k};
		end
	end
end

