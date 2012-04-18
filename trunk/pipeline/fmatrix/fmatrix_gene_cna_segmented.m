function fmatrix = fmatrix_gene_cna(segments, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = segments.Meta.Sample.ID;
if length(samples) ~= length(unique(samples))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	segments = merge_replicates(segments);
end

gcna = gene_cna(segments, 'MinRegion', 100e3);

valid = find(~any(isnan(gcna), 2));

S = length(segments.Meta.Sample.ID);
F = length(valid);

fmatrix.Samples = segments.Meta.Sample.ID;
fmatrix.Features = cell(F, 1);
fmatrix.Data = nan(F, S);

for f = 1:length(valid)
	g = valid(f);
	
	if isnan(genes.Chromosome(g))
		pos_str = '?:?:?:?';
	else
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{genes.Chromosome(g)}, ...
			genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
	end

	fmatrix.Features{f} = sprintf('N:CNA:%s:%s', genes.Name{g}, pos_str);
	fmatrix.Data(f, :) = gcna(g, :);
end

fprintf(1, '%d copy number features.\n', F);


