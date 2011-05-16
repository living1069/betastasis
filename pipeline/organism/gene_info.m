function [] = gene_info(gene_id)

global organism;
genes = organism.Genes;

g = gene_idx(gene_id);
if isnan(g), error 'Gene name was not recognized.'; end
	
fprintf(1, 'Entrez ID: %d\n', genes.EntrezID(g));

aliases = genes.Synonyms{g};
fprintf(1, 'Aliases: ');
if ~isempty(aliases)
	for k = 1:length(aliases)
		if k == length(aliases)
			fprintf(1, '%s\n', aliases{k});
		else
			fprintf(1, '%s, ', aliases{k});
		end
	end
else
	fprintf(1, '-\n');
end

chr = genes.Chromosome(g);
if ~isnan(chr)
	fprintf(1, 'Chromosome: %s (%s)\n', organism.Chromosomes.Name{chr}, ...
		genes.Strand(g));
else
	fprintf(1, 'Chromosome: -\n');
end

pos = genes.Position(g, :);
if ~any(isnan(pos))
	fprintf(1, 'Position: [%d..%d]\n', pos(1), pos(2));
else
	fprintf(1, 'Position: -\n');
end
	
txs = genes.Transcripts(g, 1:genes.TranscriptCount(g));
for t = 1:length(txs)
	tx = txs(t);
	
	tx_exons = organism.Transcripts.Exons{tx};
	exon_pos = organism.Transcripts.ExonPos{tx};
	
	fprintf(1, '\n%s:\n', organism.Transcripts.Name{tx});
	
	cds = organism.Transcripts.CDS(tx, :);
	if ~any(isnan(cds))
		tx_seq = upper(organism.Transcripts.Sequence{tx});
		fprintf(1, '- CDS [%d..%d], %s..%s\n', cds(1), cds(2), ...
			tx_seq(cds(1):cds(1)+2), tx_seq(cds(2)-2:cds(2)));
	end
	
	for x = 1:length(tx_exons)
		ex = tx_exons(x);
		fprintf(1, '- exon %s [%d..%d]\n', organism.Exons.ID{ex}, ...
			exon_pos(x, 1), exon_pos(x, 2));
	end
end

