function [] = print_exon_junction(left_exon, right_exon)

global organism;
transcripts = organism.Transcripts;
exons = organism.Exons;

left_seq = transcripts.Sequence{exons.Transcript(left_exon)};
right_seq = transcripts.Sequence{exons.Transcript(right_exon)};

left_seq = left_seq(exons.Position(left_exon, 1):exons.Position(left_exon, 2));
right_seq = right_seq(exons.Position(right_exon, 1):...
	exons.Position(right_exon, 2));

fprintf(1, '%s|%s\n', upper(left_seq(end-49:end)), upper(right_seq(1:50)));

