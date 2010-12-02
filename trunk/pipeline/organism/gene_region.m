function seq = gene_region(id, region, window)

global organism;

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

if strcmpi(region, 'TSS')
	base = nan(length(idx), 1);
	base(valid) = organism.Genes.Position(idx(valid), 1);
end

chr = nan(length(idx), 1);
chr(valid) = organism.Genes.Chromosome(idx(valid), 1);

seq = cell(length(idx), 1);
for k = 1:length(seq)
	seq{k} = '';
	if isnan(base(k)), continue, end
	
	start_offset = max(base(k) + window(1), 1);
	end_offset = min(base(k) + window(2), organism.Chromosomes.Length(chr(k)));
	
	seq{k} = organism.Chromosomes.Sequence{chr(k)}(start_offset:end_offset);
end

if ischar(id)
	seq = seq{1};
end

