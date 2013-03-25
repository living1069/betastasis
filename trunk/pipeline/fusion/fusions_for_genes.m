function [] = fusions_for_genes(reads, genesets, flank_len, varargin)

global organism;
genes = organism.Genes;
exons = organism.Exons;

S = length(reads.url);

fprintf('Constructing a list of possible exon-exon junctions...\n');

candidates = struct;
candidates.name = cell(0, 1);
candidates.exons = zeros(0, 2);
candidates.sequence = cell(0, 1);

if iscellstr(genesets), genesets = { genesets }; end

for k = 1:length(genesets)
	gset = gene_idx(genesets{k});
	gset = gset(:)';
	
	for left_gene = gset
		for right_gene = setdiff(gset, left_gene)
			left_exons = find(exons.Gene == left_gene);
			right_exons = find(exons.Gene == right_gene);
			
			for left_exon = left_exons'
				for right_exon = right_exons'
					candidates.exons(end+1, :) = [left_exon, right_exon];
				end
			end
		end
	end
end

for k = 1:size(candidates.exons, 1)
	left_exon = candidates.exons(k, 1);
	right_exon = candidates.exons(k, 2);
	
	candidates.name{k, 1} = sprintf('%s[%s]-%s[%s]', ...
		genes.Name{exons.Gene(left_exon)}, exons.ID{left_exon}, ....
		genes.Name{exons.Gene(right_exon)}, exons.ID{right_exon});
	
	left_exon_seq = exons.Sequence{left_exon};
	right_exon_seq = exons.Sequence{right_exon};

	if length(left_exon_seq) > flank_len
		left_exon_seq = left_exon_seq(end-flank_len+1:end);
	end
	if length(right_exon_seq) > flank_len
		right_exon_seq = right_exon_seq(1:flank_len);
	end
	
	candidates.sequence{k, 1} = [left_exon_seq right_exon_seq];
end

fprintf('Number of fusion junctions to check: %d\n', size(candidates.name, 1));
candidates
candidates.sequence

alignments = bowtie2_align(reads, candidates, '--score-min L,0,0');

