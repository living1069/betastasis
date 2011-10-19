function fmatrix = fmatrix_exon_expression(exon_expr, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;
exons = organism.Exons;

sample_ids = exon_expr.Meta.Sample.ID;
if length(sample_ids) ~= length(unique(sample_ids))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	exon_expr = merge_replicates(exon_expr);
end

S = size(exon_expr.Mean, 2);
E = size(exon_expr.Mean, 1);

exon_expr = exon_expr.Mean + 1;

fmatrix.Samples = sample_ids;
fmatrix.Features = cell(E, 1);
fmatrix.Data = log2(exon_expr);

for f = 1:E
	g = exons.Gene(f);
	
	if isnan(genes.Chromosome(g))
		pos_str = '?:?:?:?';
	else
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{genes.Chromosome(g)}, ...
			genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
	end
	
	fmatrix.Features{f} = sprintf('N:EXON_EXPR:%s[%s]:%s', genes.Name{g}, ...
		exons.ID{f}, pos_str);
end

fprintf(1, '%d exon expression features.\n', E);

