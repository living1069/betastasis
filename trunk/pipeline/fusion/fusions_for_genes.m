function rearrangements = fusions_for_genes(reads, genesets, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

S = length(reads.Raw);

rearrangements = struct;
rearrangements.Fusions = cell(1, S);

fprintf(1, ['Constructing a list of all possible exon-exon junctions ' ...
            'between genes of interest...\n']);

junction_seqs = {};
exon_5p = []; exon_3p = [];

if iscellstr(genesets), genesets = { genesets }; end

for k = 1:length(genesets)
	gset = gene_idx(genesets{k});
	gset = gset(:)';
	
	for left_gene = gset
		for right_gene = setdiff(gset, left_gene)
			left_exons = find(exons.Gene == left_gene);
			right_exons = find(exons.Gene == right_gene);
					
			for left_exon = left_exons
				for right_exon = right_exons
					left_exon_seq = exons.Sequence{left_exon};
					right_exon_seq = exons.Sequence{right_exon};

					if length(left_exon_seq) > 100
						left_exon_seq = left_exon_seq(end-99:end);
					end
					if length(right_exon_seq) > 100
						right_exon_seq = right_exon_seq(1:100);
					end
					
					junction_seqs{end+1, 1} = [left_exon_seq right_exon_seq];
					exon_5p(end+1, 1) = left_exon;
					exon_3p(end+1, 1) = right_exon;
				end
			end
		end
	end
end

fprintf(1, 'Number of fusion junctions to check: %d\n', num_exon_pairs);

for s = 1:S
	unaligned_read_files{k} = ptemp;
	
	fprintf(1, 'Finding unaligned transcriptome reads using Bowtie...\n');
	[status, out] = unix(sprintf(['%s/tools/bowtie/bowtie %s ' ...
		'-p4 -v2 --suppress 5,6,7,8 %s%s %s --un %s > /dev/null'], ...
		ppath, flags, bowtie_index('transcripts'), index_suffix, ...
		len_filtered_reads_tmp, unaligned_read_files{k}));
end

for k = 1:length(seq_files)
	rearrangements.Fusions{k} = validate_junctions(unaligned_read_files{k}, ...
		index_fname);
end

if isfield(reads, 'Meta')
	rearrangements.Meta = reads.Meta;
else
	rearrangements.Meta = struct;
end

rearrangements.Meta.Type = 'Genetic rearrangements';
rearrangements.Meta.Organism = organism.Name;
	







function validated = validate_junctions(unaligned_reads, index_fname)

global organism;
exons = organism.Exons;

[color, ~] = seq_read_type(unaligned_reads);

fprintf(1, 'Aligning full reads against candidate junctions...\n');
[alignments_tmp, out] = bowtie_align(unaligned_reads, index_fname, ...
	'-p4 -m4 -v2 --suppress 2,6,7,8');

fprintf(1, '%s', out);

fprintf(1, 'Reading aligned reads into memory...\n');
alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%s %s %d %s', 'Delimiter', '\t');
fclose(alignments_file);

read_ids = data{1};
junction_exons = data{2};
read_offsets = data{3} + 1;
read_sequences = data{4};
clear data;

% Calculate the lengths of the reference sequences on the 5' side of junctions.
left_lengths = zeros(length(junction_exons), 1);
for k = 1:length(junction_exons)
	ex = sscanf(junction_exons{k}, '%d,%d');
	len = exons.Position(ex(1), 2) - exons.Position(ex(1), 1) + 1;
	if len > 100, len = 100; end
	left_lengths(k) = len;
end

% Only reads that have at least 10 nucleotides on both sides of the junction
% count towards validation.
valid = false(length(read_ids), 1);
for k = 1:length(valid)
	len = length(read_sequences{k});
	offset = read_offsets(k);
	if offset <= left_lengths(k) - 7 && offset + len > left_lengths(k) + 8
		valid(k) = true;
	end
end

read_ids = read_ids(valid);
junction_exons = junction_exons(valid);
read_offsets = read_offsets(valid);
read_sequences = read_sequences(valid);

unique_junctions = unique(junction_exons);

junction_reads_count = zeros(length(unique_junctions), 1);
for k = 1:length(unique_junctions)
	junction_reads_count(k) = sum(strcmp(unique_junctions{k}, junction_exons));
end

[~, sort_indices] = sort(junction_reads_count, 1, 'descend');

validated = struct('Exons', zeros(length(unique_junctions), 2), ...
	'ReadCount', junction_reads_count);

validated.ReadSequences = cell(length(unique_junctions), ...
	max(junction_reads_count));

for k = 1:length(unique_junctions)
	validated.Exons(k, :) = sscanf(unique_junctions{k}, '%d,%d');
	
	indices = find(strcmp(unique_junctions{k}, junction_exons));
	for s = 1:length(indices)
		validated.ReadSequences{k, s} = read_sequences{indices(s)};
	end
end	


