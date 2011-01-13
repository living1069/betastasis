function [] = print_fusions(fusions, gene_blacklist)

global organism;
exons = organism.Exons;

blacklist = {};
if nargin == 2
	blacklist = gene_blacklist;
end

[~, sort_indices] = sort(fusions.ReadCount, 1, 'descend');

fprintf(1, 'Potential fusion transcripts (ordered by prevalence):\n');
for k = 1:length(sort_indices)
	idx = sort_indices(k);
	
	sequences = {};
	if isfield(fusions, 'ReadSequences')
		sequences = fusions.ReadSequences(idx, 1:fusions.ReadCount(idx));
	end
	
	print_exon_pair(fusions.Exons(idx, 1), fusions.Exons(idx, 2), ...
		fusions.ReadCount(idx), sequences, exons, blacklist);
end

return;



function [] = print_exon_pair(left_exon, right_exon, read_count, sequences, ...
	exons, blacklist)

global organism;
genes = organism.Genes;
exons = organism.Exons;

left_gene = exons.Gene(left_exon);
right_gene = exons.Gene(right_exon);

for k = 1:length(blacklist)
	if regexp(genome.Name{left_gene}, blacklist{k}), return, end
	if regexp(genome.Name{right_gene}, blacklist{k}), return, end
end

fprintf(1, '- fusion of %s[%s] and %s[%s]:\n', genes.Name{left_gene}, ...
	exons.ID{left_exon}, genes.Name{right_gene}, exons.ID{right});
fprintf(1, '  * supported by %d reads:\n', read_count);

for k = 1:length(sequences)
	fprintf(1, '    * %s\n', sequences{k});
end

return;
