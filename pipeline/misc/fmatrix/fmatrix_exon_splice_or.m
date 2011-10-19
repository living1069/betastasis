function fmatrix = fmatrix_exon_splice_or(exon_expr, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;
exons = organism.Exons;

sample_ids = exon_expr.Meta.Sample.ID;
if length(sample_ids) ~= length(unique(sample_ids))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	exon_expr = merge_replicates(exon_expr);
	sample_ids = exon_expr.Meta.Sample.ID;
end



S = size(exon_expr.Mean, 2);
E = size(exon_expr.Mean, 1);

gene_expr = nan(length(genes.Name), S);
for g = 1:length(genes.Name)
	gex = find(exons.Gene == g);
	gene_expr(g, :) = nansum(exon_expr.Mean(gex, :), 1);
end

exon_expr = exon_expr.Mean + 1;
gene_expr = gene_expr + 2;


odds_ratio = exon_expr ./ gene_expr(exons.Gene, :);
odds_ratio = odds_ratio ./ (1 - odds_ratio);

fmatrix.Samples = sample_ids;
fmatrix.Features = cell(E, 1);
fmatrix.Data = log2(odds_ratio);

for f = 1:E
	g = exons.Gene(f);
	
	if isnan(genes.Chromosome(g))
		pos_str = '?:?:?:?';
	else
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{genes.Chromosome(g)}, ...
			genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
	end
	
	fmatrix.Features{f} = sprintf('N:ESPL:%s[%s]:%s', genes.Name{g}, ...
		exons.ID{f}, pos_str);
end

fprintf(1, '%d exon expression features.\n', E);

