function rearrangements = fusions_for_genes(reads, genesets, varargin)

global organism;
exons = organism.Exons;

S = length(reads.Raw);

fprintf(1, ['Constructing a list of all possible exon-exon junctions ' ...
            'between genes of interest...\n']);

candidate_exons = zeros(0, 2);

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
					candidate_exons(end+1, :) = [left_exon, right_exon];
				end
			end
		end
	end
end

candidates.Fusions = { struct( ...
	'Exons', candidate_exons, ...
	'ReadCount', ones(size(candidate_exons, 1), 1) ...
) };

fprintf(1, 'Number of fusion junctions to check: %d\n', ...
	size(candidate_exons, 1));

rearrangements = validate_rearrangements(reads, candidates, ...
	'DiscardTxomeMatches', false, 'MaxMismatches', 1);

