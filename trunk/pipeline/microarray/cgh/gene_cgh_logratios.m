
% GENE_CGH_LOGRATIOS     Calculate aCGH logratios for individual genes
%
%    GENE_LR = GENE_CGH_LOGRATIOS(TEST, REF, PROBESETS) calculates a CGH
%    logratio for each individual gene found in organism.Genes for the currently
%    selected organism. CGH probesets must be provided in the argument
%    PROBESETS, and will determine the genomic loci to which CGH probes on the
%    array are complementary to.
%
%    Logratios are only calculated for genes whose genomic location is well
%    defined (as opposed to genes for which multiple loci have been annotated).
%    The logratios are calculated by taking the median across all probes that
%    fall within the first exon or last exon of the gene, or within
%    exons/introns between the two.
%
%    If a gene's position within the genome is not well defined, or no probes
%    target the gene, a logratio of NaN will be reported for the gene.
%
%    GENE_CGH_LOGRATIOS(..., 'MinRegion', MIN_SIZE) tells the algorithm to
%    calculate the logratios within a window of at least MIN_SIZE nucleotides.
%    If a gene is shorter than MIN_SIZE, the window within which the logratio
%    is calculated will be extended at both ends until the size MIN_SIZE is
%    reached. This option is useful for calculating the CGH logratio for short
%    genes such as microRNA genes.

% Author: Matti Annala <matti.annala@tut.fi>

function gene_lr = gene_cgh_logratios(samples, refs, probesets, varargin)

global organism;
genes = organism.Genes;

if isstruct(samples), samples = samples.Mean; end
if isstruct(refs), refs = refs.Mean; end

min_region = 0;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MinRegion')
		min_region = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

if any(size(samples) ~= size(refs))
	error 'The sample and reference matrices must have equal dimensions.';
end

S = size(samples, 2);

cnv = zeros(length(probesets.ProbeCount), S);
for k = 1:length(probesets.ProbeCount)
	probes = probesets.Probes(k, 1:probesets.ProbeCount(k));
	cnv(k, :) = median(samples(probes, :) ./ refs(probes, :), 1);
end

logratios = log2(cnv);

gene_lr = nan(length(organism.Genes.Name), size(samples, 2));

progress = Progress;

for g = 1:length(organism.Genes.Name)
	if any(isnan(organism.Genes.Position(g, :))), continue, end
		
	gene_span = genes.Position(g, :);
	if gene_span(2) - gene_span(1) < min_region
		d = min_region - (gene_span(2) - gene_span(1));
		gene_span(1) = gene_span(1) - round(d / 2);
		gene_span(2) = gene_span(2) + round(d / 2);
	end
	
	% Find all probes that target the gene we're interested in.
	idx = find(probesets.Chromosome == genes.Chromosome(g) & ...
		probesets.Offset >= gene_span(1) & ...
		probesets.Offset <= gene_span(2));
	if isempty(idx), continue, end
	
	gene_lr(g, :) = median(logratios(idx, :), 1);
	progress.update(g / length(organism.Genes.Name));
end

