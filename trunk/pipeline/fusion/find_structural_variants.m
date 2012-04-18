
% Author: Matti Annala <matti.annala@tut.fi>

function variants = find_structural_variants(reads, anchor_len, read_len)

global organism;
	
max_mismatches = 2;
min_distance = 1e6;

S = length(reads.url);

for s = 1:S
	fprintf(1, 'Analyzing sample %s...\n', reads.meta.sample_id{s});
	candidates{s} = find_anchor_based_variants(filter(reads, s), read_len, ...
		anchor_len, min_distance);
end

total_candidates = unique(cat(1, candidates{:}));

fprintf(1, 'Found a total of %d potential structural variants.\n', ...
	length(total_candidates));
	
save ~/total_candidates.mat total_candidates;

candidate_reads = repmat({{}}, length(total_candidates), S);

fprintf(1, 'Aligning complete reads to candidate junctions...\n');
alignments = bowtie_align(reads, total_candidates, ...
	sprintf('-v%d -m1', max_mismatches));
	
for s = 1:S
	for al = iterate_alignments(filter(alignments, s))
		target = str2double(al.target);
		for k = 1:length(al.read)
			t = target(k); candidate_reads{t,s}{end+1,1} = al.sequence{k};
		end
	end
end

variants.sequence = total_candidates;
variants.reads = candidate_reads;

valid = ~all(cellfun(@isempty, variants.reads), 2);
variants.sequence = variants.sequence(valid);
variants.reads = variants.reads(valid, :);


	







% FIXME: For simplification, this variant ignores anchors that align to multiple
% genomic loci.

function candidates = find_anchor_based_variants(reads, read_len, ...
	anchor_len, min_distance)
	
global organism;
chromosomes = organism.Chromosomes;

gap_len = read_len - anchor_len * 2;
if rx(reads.space{1}, 'color'), gap_len = gap_len - 2; end % FIXME
	
fprintf(1, 'Splitting reads into start and end tags...\n');
split = split_reads(reads, anchor_len);
	
fprintf(1, 'Aligning split reads to the genome...\n');
alignments = bowtie_align(split, 'genome', '-v0 -m1');
al = all_alignments(alignments);

% Here we build two vectors left_reads and right_reads, so that aligned anchors
% left_reads(k) and right_reads(k) always come from the same read.
fprintf(1, 'Searching for candidate junctions...\n');
left_reads = rx(al.read, '/1$');
read_nums = str2double(regexprep(al.read, '/[12].*$', ''));

if any(isnan(read_nums))
	error 'Not all read IDs were numeric.';
end

read_num_to_idx = nan(1, max(read_nums));
read_num_to_idx(read_nums(~left_reads)) = find(~left_reads);

left_reads = find(left_reads);
right_reads = read_num_to_idx(read_nums(left_reads));

valid = ~isnan(right_reads);
left_reads = left_reads(valid);
right_reads = right_reads(valid);

chr = chromosome_sym2num(al.target);
offset = al.offset;
strand = al.strand;

valid = find(~(chr(left_reads) == chr(right_reads) & ...
	abs(offset(left_reads) - offset(right_reads)) < min_distance));
left_reads = left_reads(valid);
right_reads = right_reads(valid);
	
fprintf(1, 'Found %d interesting junctions.\n', length(left_reads));

variant_map = containers.Map;

for k = 1:length(left_reads)
	l = left_reads(k); r = right_reads(k);
	
	%if chr(l) == chr(r) && abs(offset(l) - offset(r)) < min_distance
	%	continue;
	%end
	
	extra_bases = 5;
	
	%chr_bounds = chromosomes.Length(chr([l, r]));
	
	if strand(l) == '+'
		%[chr(l), offset(l)-extra_bases, offset(l)+anchor_len+gap_len]
		left_seq = chromosomes.Sequence{chr(l)} ...
			(offset(l)-extra_bases:offset(l)+anchor_len+gap_len);
	else
		left_seq = seqrcomplement(chromosomes.Sequence{chr(l)} ...
			(offset(l)-anchor_len-gap_len:offset(l)+extra_bases));
	end
	
	if strand(r) == '+'
		right_seq = chromosomes.Sequence{chr(r)} ...
			(offset(r)-anchor_len-gap_len:offset(r)+extra_bases);
	else
		right_seq = seqrcomplement(chromosomes.Sequence{chr(r)} ...
			(offset(r)-extra_bases:offset(r)+anchor_len+gap_len));
	end
	
	for g = 0:gap_len
		junction_seq = [left_seq(1:extra_bases+anchor_len+g) ...
			right_seq(gap_len-g+1:end)];
		variant_map(junction_seq) = 1;
	end
end

candidates = variant_map.keys()';










