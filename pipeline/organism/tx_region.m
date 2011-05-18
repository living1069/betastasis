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
	for k = na', fprintf(1, 'Could not find transcript %s.\n', id{k}); end
end

valid = ~isnan(idx);
gene_idx = transcripts.Gene(idx(valid));

window = repmat(window, length(idx), 1);
base = nan(length(idx), 1);
minus = false(length(idx), 1);

if strcmpi(region, 'TSS')
	minus(valid) = genes.Strand(gene_idx) == '-';
	window(valid & minus, :) = -window(valid & minus, end:-1:1);

	for k = find(valid)'
		tx = idx(k);
		ex = transcripts.Exons{tx};
		g = transcripts.Gene(tx);
		
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

