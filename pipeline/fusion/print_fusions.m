function [] = print_fusions(fusions, gene_blacklist)

global organism;
transcriptome = organism.Transcripts;
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
		fusions.ReadCount(idx), sequences, transcriptome, exons, blacklist);
end

return;



function [] = print_exon_pair(left_exon, right_exon, read_count, sequences, ...
	transcriptome, exons, blacklist)

global organism;
genome = organism.Genes;

left_transcript = exons.Transcript(left_exon);
right_transcript = exons.Transcript(right_exon);

left_gene = transcriptome.Gene(left_transcript);
right_gene = transcriptome.Gene(right_transcript);

for k = 1:length(blacklist)
	if regexp(genome.Name{left_gene}, blacklist{k}), return, end
	if regexp(genome.Name{right_gene}, blacklist{k}), return, end
end

fprintf(1, '- fusion of exon pair (%d, %d):\n', left_exon, right_exon);
fprintf(1, '  * exon #%d is at [%d, %d] in transcript %s of gene %s\n', ...
	left_exon, exons.Position(left_exon, 1), exons.Position(left_exon, 2), ...
	transcriptome.Name{left_transcript}, genome.Name{left_gene});
fprintf(1, '  * exon #%d is at [%d, %d] in transcript %s of gene %s\n', ...
	right_exon, exons.Position(right_exon, 1), ...
	exons.Position(right_exon, 2), ...
	transcriptome.Name{right_transcript}, genome.Name{right_gene});
fprintf(1, '  * supported by %d reads:\n', read_count);

for k = 1:length(sequences)
	fprintf(1, '    * %s\n', sequences{k});
end

return;
