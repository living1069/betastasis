
% This function goes through the given set of RNA-seq reads and looks for
% fusion gene and alternative splicing events.
%
% Inputs:
%     reads - RNA-seq reads that are to be used for discovering rearrangements.
%         The reads can be given as a filename, a cell array of filenames,
%         or a pipeline-specific data structure.
%     paired_tag_len - Length of the start and end tags.
%     prior_rearrangements - Prior rearrangements which will be incorporated
%         into the tag-based rearrangement candidates before validation. This
%         argument is optional.
% 
% Outputs:
%     probesets - Data structure describing a set of probesets, each targeting
%         a specific gene of the currently selected organism.
%
% Author: Matti Annala <matti.annala@tut.fi>

function rearrangements = find_rearrangements(reads, paired_tag_len, ...
	min_read_len, prior_rearrangements)

global organism;
transcriptome = organism.Transcripts;
exons = organism.Exons;

if nargin < 3
	min_read_len = NaN;
end

seq_files = seq_resource_files(reads);

rearrangements = struct('TagFusions', { cell(1, length(seq_files)) }, ...
                        'Fusions', { cell(1, length(seq_files)) }, ...
						'TagAltSplicings', { cell(1, length(seq_files)) }, ...
						'AltSplicings', { cell(1, length(seq_files)) });

unaligned_read_files = cell(length(seq_files), 1);
for k = 1:length(seq_files)
	len_filtered_reads_tmp = ptemp();
	unaligned_read_files{k} = ptemp();
	
	% FIXME: Length filtering doesn't currently work with colorspace reads.
	if ~isnan(min_read_len)
		fprintf(1, 'Filtering out reads shorter than %d bp...\n', min_read_len);
		[status, ~] = unix(sprintf( ...
			'%s/sources/sequencing/filter_reads_by_len.py %s %d > %s', ...
			ppath, seq_files{k}, min_read_len, len_filtered_reads_tmp));
	else
		len_filtered_reads_tmp = seq_files{k};
	end

	[flags, index_suffix] = bowtie_flags_for_reads(seq_files{k});
	
	fprintf(1, 'Finding unaligned transcriptome reads using Bowtie...\n');
	[status, out] = unix(sprintf(['%s/tools/bowtie/bowtie %s ' ...
		'-p4 -v2 --suppress 5,6,7,8 %s%s %s --un %s > /dev/null'], ...
		ppath, flags, bowtie_index('transcripts'), index_suffix, ...
		len_filtered_reads_tmp, unaligned_read_files{k}));
		
	% Remove the temporary filtered FASTA file if we have created one.
	if ~isnan(min_read_len)
		system(['rm ' len_filtered_reads_tmp]);
	end

	fprintf(1, '%s', out);
	if status ~= 0, error 'Read alignment failed.'; end
	
	[rearrangements.TagFusions{k}, rearrangements.TagAltSplicings{k}] = ...
		find_tag_rearrangements(unaligned_read_files{k}, paired_tag_len, ...
		transcriptome, exons);
end

fprintf(1, 'Pooling rearrangements for validation...\n');
total_tag_fusions = pool_rearrangements(rearrangements.TagFusions);
total_tag_alt_splicings = pool_rearrangements(rearrangements.TagAltSplicings);

if nargin == 4
	total_tag_fusions = pool_rearrangements( ...
		{ total_tag_fusions, prior_rearrangements.TagFusions{:} });
	total_tag_alt_splicings = pool_rearrangements( ...
		{ total_tag_alt_splicings, prior_rearrangements.TagAltSplicings{:} });
end

for k = 1:length(seq_files)
	rearrangements.Fusions{k} = validate_junctions(unaligned_read_files{k}, ...
		total_tag_fusions, transcriptome, exons);
	rearrangements.AltSplicings{k} = validate_junctions( ...
		unaligned_read_files{k}, total_tag_alt_splicings, transcriptome, exons);
	system(['rm ' unaligned_read_files{k}]);
end

if isfield(reads, 'Meta')
	rearrangements.Meta = reads.Meta;
else
	rearrangements.Meta = struct();
end

rearrangements.Meta.Type = 'Genetic rearrangements';
rearrangements.Meta.Organism = organism.Name;
rearrangements.Meta.GenomeVersion = organism.GenomeVersion;
rearrangements.Meta.PairedTagLength = ...
	ones(length(seq_files), 1) * paired_tag_len;
	
return;









function [tag_fusions, tag_alt_splicings] = ...
	find_tag_rearrangements(unaligned_reads, paired_tag_len, ...
	transcriptome, exons)

[flags, index_suffix] = bowtie_flags_for_reads(unaligned_reads);
[color, quality] = seq_read_type(unaligned_reads);

if quality
	error 'Reads with quality information are not supported.';
end

split_reads_tmp = ptemp();
alignments_tmp = ptemp();

color_option = 'nucleotide';
if color, color_option = 'color'; end

fprintf(1, 'Splitting unaligned reads into start and end tags...\n');
[status, ~] = unix(sprintf('%s/sources/sequencing/split_reads_into_tag_pairs.py %s %d %s > %s', ...
	ppath, unaligned_reads, paired_tag_len, color_option, split_reads_tmp));

if status ~= 0, error 'Read splitting failed.'; end



fprintf(1, 'Aligning split reads into the human exome using Bowtie...\n');
[alignments_tmp, out] = bowtie_align(split_reads_tmp, 'exons', ...
	'-p4 -m10 -v0 --suppress 5,6,7,8');

fprintf(1, '%s', out);

%fprintf(1, 'Split reads stored in %s\n', split_reads_tmp);
system(['rm ' split_reads_tmp]);

	


fprintf(1, 'Reading aligned reads into memory...\n');
alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%s %*s %d %d', 'Delimiter', '\t');
fclose(alignments_file);

system(['rm ' alignments_tmp]);

read_ids = data{1};
exon_indices = data{2};
read_offsets = data{3} + 1;
clear data;





fprintf(1, 'Searching for rearrangement events...\n');
run_ends = [ find(strcmp(read_ids(1:end-1), read_ids(2:end)) == 0); ...
	length(read_ids) ];
run_lengths = diff([0; run_ends]);
run_starts = run_ends - run_lengths + 1;

read_id_to_run = containers.Map(read_ids(run_starts), ...
	num2cell(1:length(run_lengths)));

potential_fusions = 0;
potential_alt_splicings = 0;
fusion_map = containers.Map();
alt_splicing_map = containers.Map();

progress = 0;
fprintf(1, 'Progress: 00%%');

for r = 1:length(run_starts)
	id = read_ids{run_starts(r)};
		
	if strcmp('~2', id(end-1:end)) && read_id_to_run.isKey([id(1:end-2) '~1'])
		l = read_id_to_run([id(1:end-2) '~1']);
		left_run = run_starts(l):run_ends(l);
		right_run = run_starts(r):run_ends(r);
		
		left_exons = exon_indices(left_run);
		right_exons = exon_indices(right_run);

		left_transcripts = exons.Transcript(left_exons);
		right_transcripts = exons.Transcript(right_exons);
		
		left_genes = transcriptome.Gene(left_transcripts);
		right_genes = transcriptome.Gene(right_transcripts);
		
		% Our goal is to try to find the simplest hypothesis that would explain
		% the origin of each read that did not align to the transcriptome. We
		% start by constructing a list of such transcripts where the left tag
		% maps to one exon and the right tag maps to another exon in the same
		% gene.
		if ~isempty(intersect(left_genes, right_genes))
			% Okay, at this point our hypothesis is that the start and end of
			% the read come from a known transcript, and that the middle 
			% nucleotides simply had some read errors caused by the sequencing
			% process. We can validate this hypothesis by checking if the tags
			% come from the same exon or from consecutive exons.
			
			if ~isempty(intersect(left_exons, right_exons))
				continue;    % Hypothesis: Tags come from same exon.
			end
			
			consecutive_exons = 0;
			for k = 1:length(left_exons)
				for m = 1:length(right_exons)
					if exons.Transcript(left_exons(k)) == ...
						exons.Transcript(right_exons(m)) && ...
						(exons.Position(left_exons(k), 2) == ...
						exons.Position(right_exons(m), 1) - 1 || ...
						exons.Position(left_exons(k), 1) == ...
						exons.Position(right_exons(m), 2) + 1)
						% Hypothesis: Tags come from consecutive exons.
						consecutive_exons = 1;
						break;
					end
				end
				if consecutive_exons, break, end
			end
			
			if consecutive_exons, continue, end

			% We cannot form a hypothesis for this read, so we mark all
			% exon pairs that come from a common gene as potential alternative 
			% splicing events.
			for k = 1:length(left_exons)
				for m = 1:length(right_exons)
					if transcriptome.Gene(exons.Transcript(left_exons(k))) ~=...
						transcriptome.Gene(exons.Transcript(right_exons(m)))
						continue;
					end
					
					key = sprintf('%d,%d', left_exons(k), right_exons(m));
					if ~alt_splicing_map.isKey(key)
						alt_splicing_map(key) = 0;
					end
					alt_splicing_map(key) = alt_splicing_map(key) + 1;
					
					potential_alt_splicings = potential_alt_splicings + 1;
				end
			end
		else
			% If the ~1 and ~2 reads align to some exons, but no pair of ~1
			% exon and ~2 exon belongs to the same transcript, then we have
			% reason to suspect that we're dealing with a fusion gene or
			% novel transcript.
			
			for k = 1:length(left_exons)
				for m = 1:length(right_exons)
					key = sprintf('%d,%d', left_exons(k), right_exons(m));
					if ~fusion_map.isKey(key), fusion_map(key) = 0; end
					fusion_map(key) = fusion_map(key) + 1;
						
					potential_fusions = potential_fusions + 1;
				end
			end
			
%			for k = 1:length(left_exons)
%				% If we are aligning in colorspace, then colorspace reads
%				% usually contain a primer base in the beginning and the
%				% first color is trimmed out of the read. Therefore a
%				% 20-bp read essentially becomes a 19-bp read.
%				left_end_gap = exon_len(left_exons(k)) - ...
%					(read_offsets(left_run(k)) + paired_tag_len - ...
%					colorspace) + 1;
%				
%				if left_end_gap < 0, error 'Negative left gap.'; end
%				
%				for m = 1:length(right_exons)
%					right_end_gap = read_offsets(right_run(m)) - 1;
%					total_gap = left_end_gap + right_end_gap;
%					
%					% FIXME: The size of the gap should not be hardcoded.
%					if total_gap > 5 && total_gap < 15
%						%fprintf(1, ['Potential fusion of exons (%d, %d). ' ...
%						%            'Gap: %d + %d = %d.\n'], ...
%						%	left_exons(k), right_exons(m), ...
%						%	left_end_gap, right_end_gap, ...
%						%	left_end_gap + right_end_gap);
%							
%						key = sprintf('%d,%d', left_exons(k), right_exons(m));
%						if ~fusion_map.isKey(key), fusion_map(key) = 0; end
%						fusion_map(key) = fusion_map(key) + 1;
%						
%						potential_fusions = potential_fusions + 1;
%					end
%				end
%			end
		end
	end
	
	if floor(r / length(run_starts) * 100) > progress
		progress = floor(r / length(run_starts) * 100);
		fprintf(1, '\b\b\b%02d%%', progress);
	end
end

fprintf(1, '\n');
fprintf(1, 'Found %d potential fusion events.\n', potential_fusions);
fprintf(1, 'Found %d potential alternative splicing events.\n', ...
	potential_alt_splicings);

unique_fusions = fusion_map.keys()';

tag_fusions = struct('Exons', zeros(length(unique_fusions), 2), ...
                     'ReadCount', cell2mat(fusion_map.values)');
for k = 1:length(unique_fusions)
	tag_fusions.Exons(k, :) = sscanf(unique_fusions{k}, '%d,%d');
end

unique_alt_splicings = alt_splicing_map.keys()';

tag_alt_splicings = struct('Exons', zeros(length(unique_alt_splicings), 2), ...
                       'ReadCount', cell2mat(alt_splicing_map.values)');
for k = 1:length(unique_alt_splicings)
	tag_alt_splicings.Exons(k, :) = sscanf(unique_alt_splicings{k}, '%d,%d');
end

return;











function [] = build_junction_index(index_filename, rearrangements, color, ...
	transcriptome, exons)

fasta_fname = ptemp();
fasta_file = fopen(fasta_fname, 'W');

for k = 1:size(rearrangements.Exons, 1)
	left_exon = rearrangements.Exons(k, 1);
	right_exon = rearrangements.Exons(k, 2);
	left_transcript_seq = transcriptome.Sequence{exons.Transcript(left_exon)};
	right_transcript_seq = transcriptome.Sequence{exons.Transcript(right_exon)};
	
	left_exon_seq = left_transcript_seq( ...
		exons.Position(left_exon, 1):exons.Position(left_exon, 2));
	right_exon_seq = right_transcript_seq( ...
		exons.Position(right_exon, 1):exons.Position(right_exon, 2));

	if length(left_exon_seq) > 100
		left_exon_seq = left_exon_seq(end-99:end);
	end
	if length(right_exon_seq) > 100
		right_exon_seq = right_exon_seq(1:100);
	end
	
	fprintf(fasta_file, '>%d,%d\n', left_exon, right_exon);
	fprintf(fasta_file, '%s%s\n', left_exon_seq, right_exon_seq);
end

fclose(fasta_file);

color_option = '';
index_suffix = '';
if color
	color_option = '-C';
	index_suffix = '_colorspace';
end

[status, ~] = unix(sprintf('%s/tools/bowtie/bowtie-build %s %s %s%s', ...
	ppath, color_option, fasta_fname, index_filename, index_suffix));

if status ~= 0, error 'Bowtie index construction failed.'; end

system(['rm ' fasta_fname]);
return;








function validated = validate_junctions(unaligned_reads, tag_junctions, ...
	transcriptome, exons)

junction_index_tmp = ptemp();

[color, ~] = seq_read_type(unaligned_reads);

build_junction_index(junction_index_tmp, tag_junctions, color, ...
	transcriptome, exons);

fprintf(1, 'Aligning full reads against candidate junctions...\n');
[alignments_tmp, out] = bowtie_align(unaligned_reads, junction_index_tmp, ...
	'-p4 -k10 -v2 --suppress 6,7,8');

fprintf(1, '%s', out);

system(['rm ' junction_index_tmp '*']);

fprintf(1, 'Reading aligned reads into memory...\n');
alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%s %*s %s %d %s', 'Delimiter', '\t');
fclose(alignments_file);

%fprintf(1, 'Validation alignments stored in %s\n', alignments_tmp);
system(['rm ' alignments_tmp]);

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

return;

