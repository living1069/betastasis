function gene_cna = gene_cna(samples, refs, cgh_probesets, varargin)

global organism;
genes = organism.Genes;

cna = cna_from_cgh(samples, refs, cgh_probesets, varargin{:});

gene_cna = nan(length(organism.Genes.Name), size(samples, 2));

progress = Progress;

for g = 1:length(organism.Genes.Name)
	if any(isnan(organism.Genes.Position(g, :))), continue, end
	
	% Find all probes that target the gene we're interested in.
	idx = find(cgh_probesets.Chromosome == genes.Chromosome(g) & ...
		cgh_probesets.Offset >= genes.Position(g, 1) & ...
		cgh_probesets.Offset <= genes.Position(g, 2));
	if isempty(idx), continue, end
	if any(idx ~= (idx(1):idx(end))')
		error 'Probesets are not properly ordered.';
	end
	
	gene_cna(g, :) = median(cna(idx, :), 1);
	progress.update(g / length(organism.Genes.Name));
end

