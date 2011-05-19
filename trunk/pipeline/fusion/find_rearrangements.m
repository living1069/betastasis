
% This function goes through the given set of RNA-seq reads and looks for
% fusion gene and alternative splicing events.
%
% Inputs:
%     reads - RNA-seq reads that are to be used for discovering tag_junctions.
%         The reads can be given as a filename, a cell array of filenames,
%         or a pipeline-specific data structure.
%     paired_tag_len - Length of the start and end tags.

% Author: Matti Annala <matti.annala@tut.fi>

function [rearrangements, tag_rearrangements] = find_rearrangements(reads, ...
	paired_tag_len, varargin)
	
global organism;

prior_rearrangements = [];

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'PriorRearrangements')
		prior_rearrangements = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

S = length(reads.Raw);

rearrangements.Fusions = cell(1, S);
rearrangements.AltSplicings = cell(1, S);

tag_rearrangements.Fusions = cell(1, S);
tag_rearrangements.AltSplicings = cell(1, S);

unaligned = {};

for s = 1:S
	extracted = extract_reads(filter_query(reads, s));
	
	fprintf(1, 'Finding unaligned transcriptome reads for sample #%d...\n', s);
	[al, unaligned{s}] = align_reads(extracted, 'transcripts', ...
		'MaxMismatches', 0, 'Columns', '');
	
	[tag_rearrangements.Fusions{s}, tag_rearrangements.AltSplicings{s}] = ...
		find_tag_rearrangements(unaligned{s}, paired_tag_len);
end

fprintf(1, 'Pooling rearrangements for validation...\n');
total_tag_fusions = pool_rearrangements(tag_rearrangements.Fusions);
total_tag_alt_splicings = pool_rearrangements(tag_rearrangements.AltSplicings);

if nargin == 4
	total_tag_fusions = pool_rearrangements( ...
		{ total_tag_fusions, prior_rearrangements.Fusions{:} });
	total_tag_alt_splicings = pool_rearrangements( ...
		{ total_tag_alt_splicings, prior_rearrangements.AltSplicings{:} });
end

for s = 1:S
	rearrangements.Fusions{s} = validate_junctions( ...
		unaligned{s}, total_tag_fusions);
	rearrangements.AltSplicings{s} = validate_junctions( ...
		unaligned{s}, total_tag_alt_splicings);
end

rearrangements.Meta = reads.Meta;
rearrangements.Meta.Type = 'Genetic rearrangements';
rearrangements.Meta.Organism = organism.Name;
rearrangements.Meta.PairedTagLength = ones(S, 1) * paired_tag_len;
	









function [tag_fusions, tag_alt_splicings] = ...
	find_tag_rearrangements(unaligned, paired_tag_len)
	
global organism;
exons = organism.Exons;

fprintf(1, 'Splitting unaligned reads into start and end tags...\n');
split = split_reads(unaligned, paired_tag_len);
	
fprintf(1, 'Aligning split reads to the exome using Bowtie...\n');
al = align_reads(split, 'exons', 'AllowAlignments', 3, ...
	'ReportAlignments', 3, 'MaxMismatches', 0, ...
	'Columns', 'read,strand,target,offset');

read_ids = al.ReadID;
exon_indices = str2double(al.Target);

fprintf(1, 'Searching for rearrangement events...\n');
run_ends = [ find(strcmp(read_ids(1:end-1), read_ids(2:end)) == 0); ...
	length(read_ids) ];
run_lengths = diff([0; run_ends]);
run_starts = run_ends - run_lengths + 1;

read_id_to_run = containers.Map(read_ids(run_starts), ...
	num2cell(1:length(run_lengths)));

potential_fusions = 0;
putative_alt_splicings = 0;
fusion_map = containers.Map;
alt_splicing_map = containers.Map;

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
		
		left_strands = al.Strand(left_run);
		right_strands = al.Strand(right_run);

		left_genes = exons.Gene(left_exons);
		right_genes = exons.Gene(right_exons);
		
		% Our goal is to try to find the simplest hypothesis that can explain
		% the origin of each read that did not align to the transcriptome.
		
		if ~isempty(intersect(left_genes, right_genes))
			
			% We have an alignment pair that indicates that the tags come from
			% the same gene. Hence we need not hypothesize a fusion event.
			
			if ~isempty(intersect(left_exons, right_exons))
				continue;    % Hypothesis: Tags come from the same exon.
			end
			
			%consecutive_exons = 0;
			%for k = 1:length(left_exons)
			%	for m = 1:length(right_exons)
			%		if exons.Transcript(left_exons(k)) == ...
			%			exons.Transcript(right_exons(m)) && ...
			%			(exons.Position(left_exons(k), 2) == ...
			%			exons.Position(right_exons(m), 1) - 1 || ...
			%			exons.Position(left_exons(k), 1) == ...
			%			exons.Position(right_exons(m), 2) + 1)
			%			% Hypothesis: Tags come from consecutive exons.
			%			consecutive_exons = 1;
			%			break;
			%		end
			%	end
			%	if consecutive_exons, break, end
			%end
			
			%if consecutive_exons, continue, end

			% We cannot form a hypothesis for this read, so we mark all
			% exon pairs that come from a common gene as potential alternative 
			% splicing events. That is, we assume that alternative splicing
			% events are more common than fusion events.
			for k = 1:length(left_exons)
				for m = 1:length(right_exons)
					if exons.Gene(left_exons(k)) ~= exons.Gene(right_exons(m))
						continue;
					end
					
					if left_strands(k) ~= right_strands(m), continue, end
					
					if left_strands(k) == '+'
						key = sprintf('%d,%d', left_exons(k), right_exons(m));
					elseif left_strands(k) == '-'
						key = sprintf('%d,%d', right_exons(m), left_exons(k));
					end
					
					if ~alt_splicing_map.isKey(key)
						alt_splicing_map(key) = 0;
					end
					alt_splicing_map(key) = alt_splicing_map(key) + 1;
					
					putative_alt_splicings = putative_alt_splicings + 1;
				end
			end
		else
			% If the ~1 and ~2 reads align to some exons, but no pair of ~1
			% exon and ~2 exon belongs to the same gene, then we have
			% reason to suspect that we're dealing with a fusion gene.
			
			for k = 1:length(left_exons)
				for m = 1:length(right_exons)
					
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

fprintf(1, '\n');
fprintf(1, 'Found %d potential fusion events.\n', potential_fusions);
fprintf(1, 'Found %d potential alternative splicing events.\n', ...
	putative_alt_splicings);

unique_fusions = fusion_map.keys()';
tag_fusions.Exons = zeros(length(unique_fusions), 2);
tag_fusions.ReadCount = cell2mat(fusion_map.values)';
for k = 1:length(unique_fusions)
	tag_fusions.Exons(k, :) = sscanf(unique_fusions{k}, '%d,%d');
end

unique_alt_splicings = alt_splicing_map.keys()';
tag_alt_splicings.Exons = zeros(length(unique_alt_splicings), 2);
tag_alt_splicings.ReadCount = cell2mat(alt_splicing_map.values)';
for k = 1:length(unique_alt_splicings)
	tag_alt_splicings.Exons(k, :) = sscanf(unique_alt_splicings{k}, '%d,%d');
end














function validated = validate_junctions(unaligned, tag_junctions)
	
global organism;
exons = organism.Exons;

max_read_len = 100;

junctions = struct;
junctions.Name = {};
junctions.Sequence = {};

for k = 1:size(tag_junctions.Exons, 1)
	left_exon = tag_junctions.Exons(k, 1);
	right_exon = tag_junctions.Exons(k, 2);
	
	left_exon_seq = exons.Sequence{left_exon};
	right_exon_seq = exons.Sequence{right_exon};

	if length(left_exon_seq) > max_read_len
		left_exon_seq = left_exon_seq(end-max_read_len+1:end);
	end
	if length(right_exon_seq) > max_read_len
		right_exon_seq = right_exon_seq(1:max_read_len);
	end
	
	junctions.Name{end+1, 1} = sprintf('%d,%d', left_exon, right_exon);
	junctions.Sequence{end+1, 1} = [left_exon_seq right_exon_seq];
end

fprintf(1, 'Aligning full reads against candidate junctions...\n');
al = align_reads(unaligned, junctions, 'AllowAlignments', 3, ...
	'ReportAlignments', 3, 'MaxMismatches', 0, ...
	'Columns', 'read,target,offset,sequence');

read_ids = al.ReadID;
junction_exons = al.Target;
read_sequences = al.Sequence;

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
	junction_reads_count(k) = sum(strcmp(unique_junctions{k}, junction_exons));
end

[~, sort_indices] = sort(junction_reads_count, 1, 'descend');

validated = struct;
validated.Exons = zeros(length(unique_junctions), 2);
validated.ReadCount = junction_reads_count;
validated.JunctionOffsets = zeros(length(unique_junctions), ...
	max(junction_reads_count));
validated.ReadSequences = cell(length(unique_junctions), ...
	max(junction_reads_count));

for k = 1:length(unique_junctions)
	validated.Exons(k, :) = sscanf(unique_junctions{k}, '%d,%d');
	
	indices = find(strcmp(unique_junctions{k}, junction_exons));
	for s = 1:length(indices)
		validated.JunctionOffsets(k, s) = breakpoint_offsets(indices(s));
		validated.ReadSequences{k, s} = read_sequences{indices(s)};
	end
end	

