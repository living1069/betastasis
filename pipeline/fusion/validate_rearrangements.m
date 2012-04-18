
% Author: Matti Annala <matti.annala@tut.fi>

function rearrangements = validate_rearrangements(reads, ...
	candidates, varargin)
	
global organism;
exons = organism.Exons;

max_mismatches = 0;
max_read_len = 100;

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

total_tag_fusions = pool_rearrangements(candidates.fusions);

% Construct a Bowtie index containing all of the fusion candidates.
junctions = struct;
junctions.name = {};
junctions.sequence = {};

for k = 1:size(total_tag_fusions.Exons, 1)
	left_exon = total_tag_fusions.Exons(k, 1);
	right_exon = total_tag_fusions.Exons(k, 2);
	
	left_exon_seq = exons.Sequence{left_exon};
	right_exon_seq = exons.Sequence{right_exon};

	if length(left_exon_seq) > max_read_len
		left_exon_seq = left_exon_seq(end-max_read_len+1:end);
	end
	if length(right_exon_seq) > max_read_len
		right_exon_seq = right_exon_seq(1:max_read_len);
	end
	
	junctions.name{end+1, 1} = sprintf('%d,%d', left_exon, right_exon);
	junctions.sequence{end+1, 1} = [left_exon_seq right_exon_seq];
end

tmp = temporary('fusion_validate');
fasta_tmp = [tmp 'index.fa'];
index_tmp = [tmp 'index'];

write_seq_fasta(junctions, fasta_tmp);

if length(unique(reads.space)) > 1
	error 'All reads must be in same space.';
end

color_option = '';
index_suffix = '';
if rx(reads.space{1}, 'color')
	color_option = '-C';
	index_suffix = '_colorspace';
end

index = [index_tmp index_suffix];

[status, out] = unix(sprintf('bowtie-build %s %s %s', ...
	color_option, fasta_tmp, index));
if status ~= 0, error('Bowtie index construction failed:\n%s\n.', out); end


for s = 1:S
	fprintf(['Aligning full reads against candidate junctions in ' ...
		'sample %s...\n'], reads.meta.sample_id{s});
	rearrangements.fusions{s} = validate_junctions( ...
		filter(reads, s), index_tmp, max_mismatches, max_read_len);
end

rearrangements.meta = reads.meta;
rearrangements.meta.type = 'Fusion genes';
	





function validated = validate_junctions(unaligned, index, max_mismatches, ...
	max_read_len)
	
global organism;
exons = organism.Exons;

alignments = bowtie_align(unaligned, index, ...
	sprintf('-v%d -k3 -m3', max_mismatches));

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

