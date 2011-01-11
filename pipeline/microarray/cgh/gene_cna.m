
% GENE_CNA  Calculate per-gene CNA values based on segmented copy number data
%
%    GENE_CNA = GENE_CNA(SEGMENTS, PROBESETS) calculates summarized copy number
%    alteration values GCNA for each annotated gene of the currently selected
%    organism. The CNA values are calculated based on the segmented copy number
%    data provided in argument SEGMENTS, and the CGH probesets provided in
%    PROBESETS.

% Author: Matti Annala <matti.annala@tut.fi>

function gene_cna = gene_cna(segments, cgh_probesets)

global organism;
genes = organism.Genes;

cna = cn_seg_expand(segments, cgh_probesets);

gene_cna = nan(length(organism.Genes.Name), size(cna, 2));

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

