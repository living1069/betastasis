
% GENE_CNA  Calculate per-gene CNA values based on segmented copy number data
%
%    GENE_CNA = GENE_CNA(SEGMENTS, PROBESETS) calculates summarized copy number
%    alteration values GCNA for each annotated gene of the currently selected
%    organism. The CNA values are calculated based on the segmented copy number
%    data provided in argument SEGMENTS, and the CGH probesets provided in
%    PROBESETS.

% Author: Matti Annala <matti.annala@tut.fi>

function gene_cna = gene_cna(segments, varargin)

global organism;
genes = organism.Genes;

probesets = [];
min_region = 0;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Probesets')
		probesets = varargin{k+1}; continue;
	end
	
	if strcmpi(varargin{k}, 'MinRegion')
		min_region = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isempty(probesets) && isfield(segments, 'Meta')
	probesets = platform(segments.Meta.Platform{1}, 'cgh_probesets');
end

cna = cn_seg_expand(segments, probesets);

gene_cna = nan(length(organism.Genes.Name), size(cna, 2));

progress = Progress;

for g = 1:length(genes.Name)
	if any(isnan(genes.Position(g, :))), continue, end
		
	gene_span = genes.Position(g, :);
	if gene_span(2) - gene_span(1) < min_region
		d = min_region - (gene_span(2) - gene_span(1));
		gene_span(1) = gene_span(1) - round(d / 2);
		gene_span(2) = gene_span(2) + round(d / 2);
	end
	
	% Find all probes that target the gene we're interested in.
	idx = find(probesets.Chromosome == genes.Chromosome(g) & ...
		probesets.Offset >= gene_span(1) & probesets.Offset <= gene_span(2));
	if isempty(idx), continue, end
	if any(idx ~= (idx(1):idx(end))')
		error 'Probesets are not properly ordered.';
	end
	
	gene_cna(g, :) = median(cna(idx, :), 1);
	progress.update(g / length(genes.Name));
end

