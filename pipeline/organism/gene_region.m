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
	for k = 1:length(na)
		fprintf(1, 'Could not find gene %s.\n', id{na(k)});
	end
end

valid = ~isnan(idx);

window = repmat(window, length(idx), 1);

if strcmpi(region, 'TSS')
	base = nan(length(idx), 1);
	base(valid) = genes.Position(idx(valid), 1);
	for k = 1:length(idx)
		if genes.Strand(k) == '-', window(k, :) = -window(k, end:-1:1); end
	end
end

chr = nan(length(idx), 1);
chr(valid) = genes.Chromosome(idx(valid), 1);

seq = cell(length(idx), 1);
for k = 1:length(seq)
	seq{k} = '';
	if isnan(base(k)), continue, end
		
	start_offset = max(base(k) + window(k, 1), 1);
	end_offset = min(base(k) + window(k, 2), chromosomes.Length(chr(k)));
	
	seq{k} = chromosomes.Sequence{chr(k)}(start_offset:end_offset);
end

if ischar(id)
	seq = seq{1};
end

