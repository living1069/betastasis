function pri_mirnas = find_pri_mirna()

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;
genes = organism.Genes;
exons = organism.Exons;

P = length(pre_mirnas.Name);

pri_mirnas.Type = repmat({'-'}, P, 1);
pri_mirnas.OverlapsGene = repmat({'-'}, P, 1);
pri_mirnas.Sequence = cell(P, 1);
pri_mirnas.Chromosome = nan(P, 1);
pri_mirnas.Strand = repmat(' ', P, 1);
pri_mirnas.Position = nan(P, 2);

% Recalculate gene positions based on exon coordinates. Transcript coordinates
% sometimes do not match 100% with exon coordinates.
gene_pos = nan(length(organism.Genes.Name), 2);
for g = 1:size(gene_pos, 1)
	exon_pos = exons.Position(exons.Gene == g, :);
	if isempty(exon_pos), continue, end
	gene_pos(g, :) = [min(exon_pos(:, 1)), max(exon_pos(:, 2))];
end

for p = 1:length(pre_mirnas.Name)
	% Figure out if the pre-miRNA is situated within some known gene.
	chr = pre_mirnas.Chromosome(p);
	strand = pre_mirnas.Strand(p);
	position = pre_mirnas.Position(p, :);
	
	% The pri-miRNA chromosome and strand always match the pre-miRNA ones.
	pri_mirnas.Chromosome(p) = chr;
	pri_mirnas.Strand(p) = strand;
	
	% Check if the pre-miRNA genomic position was even annotated in miRBase.
	if isnan(chr) || strand == ' ' || any(isnan(position))
		fprintf(1, 'No knowledge of the position of %s.\n', pre_mirnas.Name{p});
		pri_mirnas.Sequence{p} = pre_mirnas.Sequence{p};
		pri_mirnas.Position(p, :) = pre_mirnas.Position(p, :);
		pri_mirnas.Type{p} = 'Unknown';
		continue;
	end
	
	% Find genes that overlap with the pre-miRNA.
	neighbors = find(genes.Chromosome == chr & genes.Strand == strand);
	overlap = neighbors(gene_pos(neighbors, 1) <= position(1) & ...
		gene_pos(neighbors, 2) >= position(2));
	
	% Ignore NCBI genes that actually represent pre-miRNA positions.
	% For determining mirtrons, we are only interested in "real" genes.
	overlap = overlap(~rx(genes.Name(overlap), '^MIR'));
	
	if isempty(overlap)
		% No gene overlaps with the microRNA, so we just pick a predefined
		% amount of sequence surrounding the pre-miRNA.
		pri_pos = [position(1)-10e3, position(2)+10e3];
		
		% Note that we must watch out for gene boundaries. So we now look
		% for genes that overlap with the pri-mirna.
		overlap_5p = neighbors(gene_pos(neighbors, 2) >= pri_pos(1) & ...
			gene_pos(neighbors, 2) <= pri_pos(2));
		if ~isempty(overlap_5p)
			pri_pos(1) = max(gene_pos(overlap_5p, 2)) + 1;
		end
		
		overlap_3p = neighbors(gene_pos(neighbors, 1) >= pri_pos(1) & ...
			gene_pos(neighbors, 1) <= pri_pos(2));
		if ~isempty(overlap_3p)
			pri_pos(2) = min(gene_pos(overlap_3p, 1)) - 1;
		end
		
		pri_mirnas.Sequence{p} = organism.Chromosomes.Sequence{chr}( ...
			pri_pos(1):pri_pos(2));
		pri_mirnas.Position(p, :) = pri_pos;

		fprintf(1, '%s is an intergenic miRNA.\n', pre_mirnas.Name{p});
		
		pri_mirnas.Type{p} = 'Intergenic';
	
	else
		% Found one or more genes overlapping the microRNA.
		if length(overlap) > 1
			fprintf(1, 'More than one gene overlaps %s.\n', pre_mirnas.Name{p});
			genes.Name(overlap)
		end
		
		g = overlap(1);
		
		% Check if the pri-miRNA is in a different strand than the gene. If so,
		% we assume that the pri-miRNA has its own regulatory machinery.
		if genes.Strand(g) ~= pri_mirnas.Strand(p)
			fprintf(1, '%s is inside gene %s, but in the opposite strand.\n',...
				pre_mirnas.Name{p}, genes.Name{g});
				
			pri_mirnas.Position(p, :) = pre_mirnas.Position(p, :);
			pri_mirnas.Sequence{p} = pre_mirnas.Sequence{p};
			pri_mirnas.Type{p} = 'Intragenic (opposite strand)';
			pri_mirnas.OverlapsGene{p} = genes.Name{g};
			continue;
		end
		
		% We need to find out whether the pre-miRNA maps to an intron or exon,
		% so we build a map of exon/intron nucleotides (in logical vector form).
		gene_size = gene_pos(g, 2) - gene_pos(g, 1) + 1;
		
		exon_pos = exons.Position(exons.Gene == g, :);
		rel_exon_pos = exon_pos - gene_pos(g, 1) + 1;
		
		exonic_nucs = false(1, gene_size);
		for ex = 1:size(exon_pos, 1)
			if isnan(rel_exon_pos(ex, 1)), continue, end
			exonic_nucs(rel_exon_pos(ex, 1):rel_exon_pos(ex, 2)) = true;
		end
		
		rel_mirna_pos = position - gene_pos(g, 1) + 1;
		if any(exonic_nucs(rel_mirna_pos(1):rel_mirna_pos(2)))
			fprintf(1, '%s is inside gene %s, but overlaps with an exon.\n', ...
				pre_mirnas.Name{p}, genes.Name{overlap(1)});
			pri_mirnas.Position(p, :) = pre_mirnas.Position(p, :);
			pri_mirnas.Sequence{p} = pre_mirnas.Sequence{p};
			pri_mirnas.Type{p} = 'Intragenic (overlaps exon)';
			pri_mirnas.OverlapsGene{p} = genes.Name{g};
			continue;
		end
		
		% Find the 5' end of the intron.
		pri_pos = nan(1, 2);
		for k = rel_mirna_pos(1):-1:1
			if exonic_nucs(k) == true
				pri_pos(1) = (gene_pos(g, 1) - 1) + k + 1;
				break;
			end
		end
		
		% Find the 3' end of the intron.
		for k = rel_mirna_pos(2):length(exonic_nucs)
			if exonic_nucs(k) == true
				pri_pos(2) = (gene_pos(g, 1) - 1) + k - 1;
				break;
			end
		end
		
		% This pri-miRNA looks like a valid mirtron.
		fprintf(1, '%s is a mirtron transcribed from %s.\n', ...
			pre_mirnas.Name{p}, genes.Name{g});

		pri_mirnas.Sequence{p} = organism.Chromosomes.Sequence{chr}( ...
			pri_pos(1):pri_pos(2));
		pri_mirnas.Position(p, :) = pri_pos;
		pri_mirnas.Type{p} = 'Mirtron';
		pri_mirnas.OverlapsGene{p} = genes.Name{g};
	end
end

mirtron = strcmp('Mirtron', pri_mirnas.Type);
sum(mirtron)

pri_mirnas.Matures = cell(length(pri_mirnas.Type), 1);
for p = 1:length(pri_mirnas.Type)
	pri_mirnas.Matures{p} = pre_mirnas.Matures(p, 1:pre_mirnas.MatureCount(p));
end


if 0
[uniq_pri, idx, rev] = unique(pri_mirnas.Sequence);
save ~/rev.mat rev;

for p = 1:length(uniq_pri)
	pri_mirnas.Premirnas{p, 1} = find(rev == p);
end

pri_mirnas.Name = cell(length(uniq_pri), 1);
for p = 1:length(uniq_pri)
	pri_pre = pri_mirnas.Premirnas{p}(1);
	pri_mirnas.Name{p} = regexprep(pre_mirnas.Name{pri_pre}, ...
		'(mir|let)-', 'pri$1-');
end

pri_mirnas.Type = pri_mirnas.Type(idx);
pri_mirnas.Sequence = pri_mirnas.Sequence(idx);
pri_mirnas.Chromosome = pri_mirnas.Chromosome(idx);
pri_mirnas.Strand = pri_mirnas.Strand(idx);
pri_mirnas.Position = pri_mirnas.Position(idx, :);
pri_mirnas.OverlapsGene = pri_mirnas.OverlapsGene(idx);
end

