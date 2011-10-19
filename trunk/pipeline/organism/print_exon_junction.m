function [] = print_exon_junction(exon_5p, exon_3p)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

tokens = regexpi(exon_5p, '(.+)\[(.+)\]', 'tokens');
if length(tokens) ~= 1, error 'Invalid 5'' exon name.'; end
token = tokens{1}; gene_5p = token{1}; exon_5p = token{2};
	
tokens = regexpi(exon_3p, '(.+)\[(.+)\]', 'tokens');
if length(tokens) ~= 1, error 'Invalid 3'' exon name.'; end
token = tokens{1}; gene_3p = token{1}; exon_3p = token{2};

left_gene = gene_idx(gene_5p);
right_gene = gene_idx(gene_3p);

if isnan(left_gene), error 'Unknown 5'' gene.'; end
if isnan(right_gene), error 'Unknown 3'' gene.'; end

exons_5p = [];
for tx = genes.Transcripts(left_gene, 1:genes.TranscriptCount(left_gene))
	exons_5p = [exons_5p; transcripts.Exons{tx}];
end
exons_3p = [];
for tx = genes.Transcripts(right_gene, 1:genes.TranscriptCount(right_gene))
	exons_3p = [exons_3p; transcripts.Exons{tx}];
end

exons_5p = unique(exons_5p);
exons_3p = unique(exons_3p);

left_exon = exons_5p(strcmp(exon_5p, exons.ID(exons_5p)));
right_exon = exons_3p(strcmp(exon_3p, exons.ID(exons_3p)));

left_seq = exons.Sequence{left_exon};
right_seq = exons.Sequence{right_exon};

fprintf(1, '%s|%s\n', upper(left_seq), upper(right_seq));

