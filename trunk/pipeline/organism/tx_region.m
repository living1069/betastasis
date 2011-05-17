function seq = tx_region(id, region, window)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;
chromosomes = organism.Chromosomes;

if numel(window) ~= 2 || ~isnumeric(window) || window(1) > window(2)
	error 'The window must be given as a two element vector.';
end

if isnumeric(id)
	idx = id;
else
	idx = transcript_idx(id);
	na = find(isnan(idx));
	for k = 1:length(na)
		fprintf(1, 'Could not find transcript %s.\n', id{na(k)});
	end
end

valid = find(~isnan(idx));
gene_idx = transcripts.Gene(idx);

window = repmat(window, length(idx), 1);
base = nan(length(idx), 1);

if strcmpi(region, 'TSS')
	neg = genes.Strand(gene_idx) == '-';
	window(valid(neg), :) = -window(valid(neg), end:-1:1);

	for k = 1:length(idx)
		if ~valid(k), continue, end
		
		tx = idx(k);
		ex = transcripts.Exons{tx};
		g = transcripts.Gene(tx);
		
		genes.Name{g}
		
		if genes.Strand(g) == ' '
			error('Gene %s is in unknown strand.', genes.Name{g});
		elseif genes.Strand(g) == '+'
			base(k) = min(exons.Position(ex, 1));
		elseif genes.Strand(g) == '-'
			base(k) = max(exons.Position(ex, 2));
		end
	end
end

chr = nan(length(idx), 1);
chr(valid) = genes.Chromosome(gene_idx, 1);

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

