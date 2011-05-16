function [] = print_transcript_exons(transcript)

global organism;

if ischar(transcript)
	transcript = transcript_idx(transcript);
end

exons = organism.Transcripts.Exons{transcript};
exon_pos = organism.Transcripts.ExonPos{transcript};

for k = 1:length(exons)
	fprintf(1, 'Exon %s [%d..%d]:\n', organism.Exons.ID{exons(k)}, ...
		exon_pos(k, 1), exon_pos(k, 2));
	fprintf(1, '%s\n\n', organism.Exons.Sequence{exons(k)});
end

