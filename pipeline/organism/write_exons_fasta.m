function [] = write_exons_fasta(filename, options)

global organism;
transcriptome = organism.Transcripts;
exons = organism.Exons;

only_unique = 0;
if nargin == 2 && strcmp('unique', options)
	only_unique = 1;
end

seen_exons = containers.Map();

output = fopen(filename, 'w');

for k = 1:length(exons.Transcript)
	transcript_seq = transcriptome.Sequence{exons.Transcript(k)};
	exon_seq = transcript_seq(exons.Position(k, 1):exons.Position(k, 2));
	
	if only_unique
		if seen_exons.isKey(exon_seq), continue, end
		seen_exons(exon_seq) = 1;
	end
	
	fprintf(output, '>%d\n%s\n', k, exon_seq);
end

fclose(output);

