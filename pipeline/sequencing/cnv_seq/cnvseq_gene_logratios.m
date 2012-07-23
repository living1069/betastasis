
% Author: Matti Annala <matti.annala@tut.fi>

function gene_cna = cnvseq_gene_logratios(cna, varargin)

global organism;
genes = organism.Genes;

min_region = 0;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'min.*region')
		min_region = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

S = size(cna.mean, 2);
gene_cna = struct;
gene_cna.meta = cna.meta;
gene_cna.rows = struct;
gene_cna.rows.gene = genes.Name;
gene_cna.mean = nan(length(genes.Name), S);

for g = 1:length(genes.Name)
	if any(isnan(genes.Position(g, :))), continue, end
		
	gene_span = genes.Position(g, :);
	if gene_span(2) - gene_span(1) < min_region
		d = min_region - (gene_span(2) - gene_span(1));
		gene_span(1) = gene_span(1) - round(d / 2);
		gene_span(2) = gene_span(2) + round(d / 2);
	end
	
	% Find all windows that target the gene we're interested in.
	idx = find(cna.rows.chromosome == genes.Chromosome(g) & ...
		cna.rows.position >= gene_span(1) & cna.rows.position <= gene_span(2));
	if isempty(idx), continue, end
		
	gene_cna.mean(g, :) = nanmedian(cna.mean(idx, :), 1);
end

