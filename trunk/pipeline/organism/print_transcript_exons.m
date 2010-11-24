function [] = print_transcript_exons(transcript)

global organism;
seq = organism.Transcripts.Sequence{transcript};
exon_indices = find(organism.Exons.Transcript == transcript);

for k = 1:length(exon_indices)
	idx = exon_indices(k);
	fprintf(1, 'Exon %d [%d..%d]:\n', k, organism.Exons.Position(idx, 1), ...
		organism.Exons.Position(idx, 2));
	fprintf(1, '%s\n\n', seq(organism.Exons.Position(idx, 1): ...
		organism.Exons.Position(idx, 2)));
end

