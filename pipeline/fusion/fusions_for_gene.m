function rearrangements = fusions_for_gene(reads, gene, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

gene_role = '5p';
max_mismatches = 2;
max_read_len = 100;

S = length(reads.Raw);

if ischar(gene)
	gene = gene_idx(gene);
end

rearrangements = struct;
rearrangements.Fusions = cell(1, S);

fprintf(1, ['Constructing a list of all possible exon-exon junctions ' ...
    'between exons of %s and other human exons...\n'], ...
	genes.Name{gene});

gene_exons = find(organism.Exons.Gene == gene);
other_exons = find(organism.Exons.Gene ~= gene);

junctions = struct;
junctions.Exons = nan(length(gene_exons) * length(other_exons), 2);
junctions.Sequence = cell(length(gene_exons) * length(other_exons), 1);

fprintf(1, 'Number of fusion junctions to check: %d\n', ...
	length(gene_exons) * length(other_exons));

J = 0;

for gex = gene_exons'
	for oex = other_exons'
		if strcmpi(gene_role, '5p')
			left_exon_seq = exons.Sequence{gex};
			right_exon_seq = exons.Sequence{oex};
			
			if length(left_exon_seq) > max_read_len
				left_exon_seq = left_exon_seq(end-max_read_len+1:end);
			end
			if length(right_exon_seq) > max_read_len
				right_exon_seq = right_exon_seq(1:max_read_len);
			end
		end
			
		J = J + 1;
		junctions.Sequence{J} = [left_exon_seq right_exon_seq];
		junctions.Exons(J, :) = [gex oex];
	end
end

for s = 1:S
	fprintf(1, 'Finding unaligned transcriptome reads for sample #%d...\n', s);
	[~, unaligned] = align_reads(filter_query(reads, s), 'transcripts', ...
		'MaxMismatches', max_mismatches, 'Columns', '');

	fprintf(1,'Aligning full unaligned reads against candidate junctions...\n');
	al = align_reads(unaligned, junctions, ...
		'AllowAlignments', 1, 'MaxMismatches', max_mismatches, ...
		'Columns', 'read,target,offset,sequence');
		
	read_ids = al.ReadID;
	junction_exons = str2double(al.Target);
	read_sequences = al.Sequence;

	valid = false(length(read_ids), 1);
	breakpoint_offsets = zeros(length(read_ids), 1);
	for k = 1:length(junction_exons)
		% Calculate the length of the validation sequence on the 5' side of the
		% junction.
		ex = junctions.Exons(junction_exons(k), :);
		left_len = length(exons.Sequence{ex(1)});
		if left_len > max_read_len, left_len = max_read_len; end

		% Only reads that have at least 5 nucleotides on both sides of the
		% junction count towards validation.
		len = length(read_sequences{k});
		offset = al.Offset(k);
		if offset <= left_len - 5 && offset + len > left_len + 5
			valid(k) = true;
		end
		
		breakpoint_offsets(k) = left_len - offset + 1;
	end

	%read_ids = read_ids(valid);
	junction_exons = junction_exons(valid);
	%al.Offset = al.Offset(valid);
	breakpoint_offsets = breakpoint_offsets(valid);
	read_sequences = read_sequences(valid);

	unique_junctions = unique(junction_exons);

	junction_reads_count = zeros(length(unique_junctions), 1);
	for k = 1:length(unique_junctions)
		junction_reads_count(k) = sum(unique_junctions(k) == junction_exons);
	end

	[~, sort_indices] = sort(junction_reads_count, 1, 'descend');

	validated = struct;
	validated.Exons = zeros(length(unique_junctions), 2);
	validated.ReadCount = junction_reads_count;
	validated.JunctionOffsets = cell(length(unique_junctions), 1);
	validated.ReadSequences = cell(length(unique_junctions), 1);

	for k = 1:length(unique_junctions)
		validated.Exons(k, :) = junctions.Exons(unique_junctions(k), :);
		
		indices = find(unique_junctions(k) == junction_exons);
		validated.JunctionOffsets{k} = breakpoint_offsets(indices);
		validated.ReadSequences{k} = read_sequences(indices);
	end	

	rearrangements.Fusions{s} = validated;
end

rearrangements.Meta = reads.Meta;
rearrangements.Meta.Type = 'Genetic rearrangements';
rearrangements.Meta.Organism = organism.Name;

