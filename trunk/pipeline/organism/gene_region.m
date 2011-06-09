
% GENE_REGION     Retrieve genomic sequence close to a gene feature
%
%   SEQ = GENE_REGION(

function seq = gene_region(id, region, window)

global organism;
genes = organism.Genes;
chromosomes = organism.Chromosomes;

if numel(window) ~= 2 || ~isnumeric(window) || window(1) > window(2)
	error 'The window must be given as a two element vector.';
end

if isnumeric(id)
	idx = id;
else
	idx = gene_idx(id);
	na = find(isnan(idx));
	for k = na', fprintf(1, 'Could not find gene %s.\n', id{k}); end
end

valid = ~isnan(idx);

window = repmat(window, length(idx), 1);
base = nan(length(idx), 1);
minus = false(length(idx), 1);

if strcmpi(region, 'TSS') || strcmpi(region, 'start')
	minus(valid) = genes.Strand(idx(valid)) == '-';
	window(valid & minus, :) = -window(valid & minus, end:-1:1);
	base(valid) = genes.Position(idx(valid), 1);
	base(valid & minus) = genes.Position(idx(valid), 2);
	
	valid(isnan(base)) = false;
end

chr = nan(length(idx), 1);
chr(valid) = genes.Chromosome(idx(valid), 1);

seq = repmat({''}, length(idx), 1);
for k = find(valid)'
	start_offset = max(base(k) + window(k, 1), 1);
	end_offset = min(base(k) + window(k, 2), chromosomes.Length(chr(k)));
	
	seq{k} = chromosomes.Sequence{chr(k)}(start_offset:end_offset);
	if minus(k), seq{k} = seqrcomplement(seq{k}); end
end

if ischar(id)
	seq = seq{1};
end

