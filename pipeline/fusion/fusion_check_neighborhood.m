function good = fusion_check_neighborhood(fusion_names, fusions, distance)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;
exons = organism.Exons;

if nargin < 3, distance = 50e3; end

if ischar(fusion_names), fusion_names = { fusion_names }; end
	
good = true(length(fusion_names), 1);

if isfield(fusions, 'Fusions')
	fusions = fusions.Fusions;
elseif isstruct(fusions)
	fusions = { fusions };
end

S = length(fusions);

chr_seq = cell(length(chromosomes.Sequence), 1);

for f = 1:length(fusion_names)
	
	tokens = regexpi(fusion_names{f}, '(.+):(.+)', 'tokens');
	if length(tokens) ~= 1, error 'Invalid fusion specifier.'; end

	token = tokens{1};
	gene_5p_str = token{1}; gene_3p_str = token{2};

	gene_5p = gene_idx(gene_5p_str);
	gene_3p = gene_idx(gene_3p_str);

	if isnan(gene_5p), error('Unknown 5'' gene %s.', gene_5p_str); end
	if isnan(gene_3p), error('Unknown 3'' gene %s.', gene_3p_str); end
		
	chr_5p = genes.Chromosome(gene_5p);
	chr_3p = genes.Chromosome(gene_3p);
	range_5p = genes.Position(gene_5p, :);
	range_3p = genes.Position(gene_3p, :);
	
	if isnan(chr_5p) || any(isnan(range_5p)) || isnan(chr_3p) || ...
		any(isnan(range_3p))
		fprintf(1, 'Gene loci not known for %s. Skipping...\n', ...
			fusion_names{f});
		continue;
	end
	
	% If the genes are in the same chromosome, we have to check if the genes
	% are really close to one another.
	if chr_5p == chr_3p && ...
		(range_5p(1) <= range_3p(1) && range_5p(2) >= range_3p(1)) || ...
		(range_5p(1) >= range_3p(2) && range_5p(2) <= range_3p(2))
		fprintf(1, ['Genes for %s are in too close vicinity. ' ...
			'Cannot check for sequence homology.\n'], fusion_names{f});
	end
	
	range_5p = [range_5p(1)-distance, range_5p(2)+distance];
	range_3p = [range_3p(1)-distance, range_3p(2)+distance];
	
	range_5p = max(range_5p, 1);
	range_5p = min(range_5p, chromosomes.Length(chr_5p));
	
	range_3p = max(range_3p, 1);
	range_3p = min(range_3p, chromosomes.Length(chr_3p));
	
	locus_seq_5p = chromosomes.Sequence{chr_5p}(range_5p(1):range_5p(2));
	locus_seq_3p = chromosomes.Sequence{chr_3p}(range_3p(1):range_3p(2));
	
	% Now we try to find the longest anchors on both sides of the junction,
	% for each different exon pair.
	fusion_map = containers.Map;
	flanks_5p = {};
	flanks_3p = {};
	
	for s = 1:S
		for e = 1:size(fusions{s}.Exons, 1)
			left_exon = fusions{s}.Exons(e, 1);
			right_exon = fusions{s}.Exons(e, 2);
			
			left_gene = exons.Gene(left_exon);
			right_gene = exons.Gene(right_exon);
			
			if ~(left_gene == gene_5p && right_gene == gene_3p), continue, end
				
			fusion_id = sprintf('%d,%d', left_exon, right_exon);
			if ~fusion_map.isKey(fusion_id)
				flanks_5p{end+1, 1} = repmat(' ', 0, 100);
				flanks_3p{end+1, 1} = repmat(' ', 0, 100);
				fusion_map(fusion_id) = length(flanks_5p);
			end
			idx = fusion_map(fusion_id);
			
			seqs = fusions{s}.ReadSequences{e};
			junc_offsets = fusions{s}.JunctionOffsets{e};
			
			for k = 1:length(seqs)
				left_len = junc_offsets(k);
				right_len = length(seqs{k}) - left_len;
				
				flanks_5p{idx}(end+1, end-left_len+1:end) = seqs{k}(1:left_len);
				flanks_3p{idx}(end+1, 1:right_len) = seqs{k}(left_len+1:end);
			end
		end
	end
	
	nucs = 'ACGT';
	
	for e = 1:length(flanks_5p)
		% Calculate the consensus 5p flank.
		flanks = flanks_5p{e};
		acgt = zeros(4, size(flanks, 2));
		acgt(1, :) = sum(flanks == 'A', 1);
		acgt(2, :) = sum(flanks == 'C', 1);
		acgt(3, :) = sum(flanks == 'G', 1);
		acgt(4, :) = sum(flanks == 'T', 1);
		acgt = acgt(:, sum(acgt, 1) >= 1);
		
		[~, idx] = max(acgt, [], 1);
		consensus_5p = nucs(idx);
		
		% Calculate the consensus 3p flank.
		flanks = flanks_3p{e};
		acgt = zeros(4, size(flanks, 2));
		acgt(1, :) = sum(flanks == 'A', 1);
		acgt(2, :) = sum(flanks == 'C', 1);
		acgt(3, :) = sum(flanks == 'G', 1);
		acgt(4, :) = sum(flanks == 'T', 1);
		acgt = acgt(:, sum(acgt, 1) >= 1);
	
		[~, idx] = max(acgt, [], 1);
		consensus_3p = nucs(idx);
		
		%consensus_5p
		%consensus_3p
		
		% Now we need to check if we can find the consensus flank in the genomic
		% neighborhood of the other flank.
		pos = strfind(locus_seq_5p, consensus_3p);
		pos = [pos, strfind(locus_seq_5p, seqrcomplement(consensus_3p))];
		if ~isempty(pos)
			fprintf(1, 'Found the 3'' anchor close to 5'' gene.\n');
			good(f) = false;
		end
		
		pos = strfind(locus_seq_3p, consensus_5p);
		pos = [pos, strfind(locus_seq_3p, seqrcomplement(consensus_5p))];
		if ~isempty(pos)
			fprintf(1, 'Found the 5'' anchor close to 3'' gene.\n');
			good(f) = false;
		end
	end
end




