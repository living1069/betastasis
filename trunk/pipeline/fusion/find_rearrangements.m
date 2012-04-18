
% This function goes through the given set of RNA-seq reads and looks for
% fusion gene and alternative splicing events.
%
% Inputs:
%     reads - RNA-seq reads that are to be used for discovering tag_junctions.
%         The reads can be given as a filename, a cell array of filenames,
%         or a pipeline-specific data structure.
%     paired_tag_len - Length of the start and end tags.

% Author: Matti Annala <matti.annala@tut.fi>

function rearrangements = find_rearrangements(txome_unaligned_reads, ...
	paired_tag_len, varargin)
	
reads = txome_unaligned_reads;
	
prior_fusions = [];
max_mismatches = 2;

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if rx(varargin{k}, 'max.*mismatch')
		max_mismatches = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

S = length(reads.url);

rearrangements.fusions = cell(1, S);

for s = 1:S
	fprintf('Searching for candidate fusions in sample %s...\n', ...
		reads.meta.sample_id{s});
	tag_fusions = find_tag_rearrangements(filter(reads, s), paired_tag_len);
		
	fprintf(['Aligning full reads against candidate junctions ' ...
		'in sample %s...\n'], reads.meta.sample_id{s});
	rearrangements.fusions{s} = validate_junctions( ...
		filter(reads, s), tag_fusions, max_mismatches);
end

rearrangements.meta = reads.meta;
rearrangements.meta.type = 'Fusion genes';
rearrangements.paired_tag_length = ones(1, S) * paired_tag_len;
	









function tag_fusions = find_tag_rearrangements(unaligned, paired_tag_len)
	
global organism;
exons = organism.Exons;

fprintf('Splitting unaligned reads into start and end tags...\n');
split = split_reads(unaligned, paired_tag_len);
	
fprintf('Aligning split reads to the exome using Bowtie...\n');
alignments = bowtie_align(split, 'exons', '-k3 -m3 -v0');

al = all_alignments(alignments);
read_ids = al.read;
exon_indices = str2double(al.target);

fprintf('Searching for rearrangement events...\n');
run_ends = [ find(strcmp(read_ids(1:end-1), read_ids(2:end)) == 0); ...
	length(read_ids) ];
run_lengths = diff([0; run_ends]);
run_starts = run_ends - run_lengths + 1;

read_id_to_run = containers.Map(read_ids(run_starts), ...
	num2cell(1:length(run_lengths)));

potential_fusions = 0;
fusion_map = containers.Map;

progress = Progress;

for r = 1:length(run_starts)
	id = read_ids{run_starts(r)};
	
	if ~strcmp('/2', id(end-1:end)), continue, end
	
	if read_id_to_run.isKey([id(1:end-2) '/1'])
		l = read_id_to_run([id(1:end-2) '/1']);
		left_run = run_starts(l):run_ends(l);
		right_run = run_starts(r):run_ends(r);
		
		left_exons = exon_indices(left_run);
		right_exons = exon_indices(right_run);
		
		left_strands = al.strand(left_run);
		right_strands = al.strand(right_run);

		left_genes = exons.Gene(left_exons);
		right_genes = exons.Gene(right_exons);
		
		% Our goal is to try to find the simplest hypothesis that can explain
		% the origin of each read that did not align to the transcriptome.
		
		if isempty(intersect(left_genes, right_genes))
			% If the ~1 and ~2 reads align to some exons, but no pair of ~1
			% exon and ~2 exon belongs to the same gene, then we have
			% reason to suspect that we're dealing with a fusion gene.
			
			for k = 1:length(left_exons)
				%if blacklist(left_genes(k)), continue, end
					
				for m = 1:length(right_exons)
					%if blacklist(right_genes(m)), continue, end
					
					% FIXME: Add support for inversions.
					if left_strands(k) ~= right_strands(m), continue, end
						
					if left_strands(k) == '+'
						key = sprintf('%d,%d', left_exons(k), right_exons(m));
					elseif left_strands(k) == '-'
						key = sprintf('%d,%d', right_exons(m), left_exons(k));
					end
					
					if ~fusion_map.isKey(key), fusion_map(key) = 0; end
					fusion_map(key) = fusion_map(key) + 1;
					
					potential_fusions = potential_fusions + 1;
				end
			end
		end
	end
	
	progress.update(r / length(run_starts));
end

fprintf('\n');
fprintf('Found %d potential fusion events.\n', potential_fusions);

unique_fusions = fusion_map.keys()';
tag_fusions.Exons = zeros(length(unique_fusions), 2);
tag_fusions.ReadCount = cell2mat(fusion_map.values)';
for k = 1:length(unique_fusions)
	tag_fusions.Exons(k, :) = sscanf(unique_fusions{k}, '%d,%d');
end















function validated = validate_junctions(unaligned, tag_fusions, max_mismatches)
	
global organism;
exons = organism.Exons;

max_read_len = 100;

junctions = struct;
junctions.name = cell(size(tag_fusions.Exons, 1), 1);
junctions.sequence = cell(size(tag_fusions.Exons, 1), 1);

for k = 1:size(tag_fusions.Exons, 1)
	left_exon = tag_fusions.Exons(k, 1);
	right_exon = tag_fusions.Exons(k, 2);
	
	left_exon_seq = exons.Sequence{left_exon};
	right_exon_seq = exons.Sequence{right_exon};

	if length(left_exon_seq) > max_read_len
		left_exon_seq = left_exon_seq(end-max_read_len+1:end);
	end
	if length(right_exon_seq) > max_read_len
		right_exon_seq = right_exon_seq(1:max_read_len);
	end
	
	junctions.name{k} = sprintf('%d,%d', left_exon, right_exon);
	junctions.sequence{k} = [left_exon_seq right_exon_seq];
end

alignments = bowtie_align(unaligned, junctions, ...
	sprintf('-k3 -m3 -v%d', max_mismatches));
	
al = all_alignments(alignments);
read_ids = al.read;
junction_exons = al.target;
read_sequences = al.sequence;

valid = false(length(read_ids), 1);
breakpoint_offsets = zeros(length(read_ids), 1);
for k = 1:length(junction_exons)
	% Calculate the length of the validation sequence on the 5' side of the
	% junction.
	ex = sscanf(junction_exons{k}, '%d,%d');
	left_len = length(exons.Sequence{ex(1)});
	if left_len > max_read_len, left_len = max_read_len; end

	% Only reads that have at least 5 nucleotides on both sides of the junction
	% count towards validation.
	len = length(read_sequences{k});
	offset = al.offset(k);
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
	junction_reads_count(k) = sum(strcmp(unique_junctions{k}, junction_exons));
end

[~, sort_indices] = sort(junction_reads_count, 1, 'descend');

validated = struct;
validated.Exons = zeros(length(unique_junctions), 2);
validated.ReadCount = junction_reads_count;
validated.JunctionOffsets = cell(length(unique_junctions), 1);
validated.ReadSequences = cell(length(unique_junctions), 1);

for k = 1:length(unique_junctions)
	validated.Exons(k, :) = sscanf(unique_junctions{k}, '%d,%d');
	
	indices = find(strcmp(unique_junctions{k}, junction_exons));
	validated.JunctionOffsets{k} = breakpoint_offsets(indices);
	validated.ReadSequences{k} = read_sequences(indices);
end	

