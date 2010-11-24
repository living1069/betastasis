function [] = exome_stats()

global organism;
exons = organism.Exome;

exon_seq = cell(length(exons.Transcript), 1);
for k = 1:length(exon_seq)
	s = organism.Transcripts.Sequence{exons.Transcript(k)};
	exon_seq{k} = s(exons.Position(k, 1):exons.Position(k, 2));
end

fprintf(1, 'Total number of exons: %d\n', length(exons.Transcript));
fprintf(1, 'Unique exons: %d\n', length(unique(exon_seq)));

exon_count = zeros(length(organism.Transcripts.Name), 1);
for k = 1:length(exon_count)
	exon_count(k) = sum(exons.Transcript == k);
end

fprintf(1, 'Exons per transcript (mean): %f\n', mean(exon_count));
fprintf(1, 'Exons per transcript (median): %f\n', median(exon_count));

