
% READS_TO_IGV      Construct IGV tracks of genomic read alignments
%
%    READS_TO_IGV(READS, FILE_PREFIX) aligns the short sequence reads in dataset
%    READS against the genome, and constructs IGV tracks for each sample. The
%    IGV tracks will be stored with filenames based on FILE_PREFIX, so that
%    if FILE_PREFIX = 'read_track', the tracks will be named 'read_track_1.igv',
%    'read_track_2.igv' and so forth.
%
%    READS_TO_IGV(..., 'Granularity', GR) sets the size of the bins to which
%    the genome will be divided up when performing run length encoding to
%    decrease the size of the IGV tracks. Default is 10 bases.
%    
%    Alignments against the genome are done allowing a maximum of 2 base
%    mismatches in the reads. A maximum of 10 alignments are reported for each
%    read, and reads with more than 10 alignments are not reported.
%    See HELP ALIGN_READS for the list of additional alignment options.

% Author: Matti Annala <matti.annala@tut.fi>

function [] = reads_to_igv(reads, igv_file_prefix, varargin)

global organism;

granularity = 10;

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Granularity')
		granularity = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

seq_files = seq_resource_files(reads);

track_map = containers.Map(strcat('chr', organism.Chromosomes.Name)', ...
	num2cell(1:length(organism.Chromosomes.Name)));

for seq_file = 1:length(seq_files)
	fprintf(1, 'Aligning reads to genome...\n');
	al = align_reads(seq_files{seq_file}, 'genome', ...
		'Columns', 'target,offset,sequence', ...
		'MaxMismatches', 2, 'AllowAlignments', 10, varargin{:});

	chromosomes = al.Target;
	start_pos = int32(al.Offset);
	
	seqlens = zeros(length(al.Sequence), 1, 'int32');
	for k = 1:length(seqlens)
		seqlens(k) = length(al.Sequence{k});
	end

	track_indices = cell2mat(track_map.values(chromosomes));
	start_slots = idivide(start_pos, granularity) + 1;
	end_slots = idivide(start_pos + seqlens - 1, granularity) + 1;

	tracks = zeros(25, max(end_slots), 'single');

	fprintf(1, 'Constructing chromosome score tracks...\n');

	for k = 1:length(start_slots)
		slot_range = start_slots(k):end_slots(k);
		tracks(track_indices(k), slot_range) = ...
			tracks(track_indices(k), slot_range) + ones(1, length(slot_range));
	end

	igv_file = sprintf('%s_%d.igv', igv_file_prefix, seq_file);
	out = fopen(igv_file, 'W');
	fprintf(out, 'Chromosome\tStart\tEnd\tFeature\tReads\n');

	for t = 1:25
		chromosome = organism.Chromosomes.Name{t};

		% This does run-length encoding (RLE) in order to optimize the .igv file
		run_ends = [ find(tracks(t, 1:end-1) ~= tracks(t, 2:end)) ...
			size(tracks, 2) ];
		run_lengths = diff([0 run_ends]);
		run_scores = tracks(t, run_ends);
		
		pos = 0;
		for r = 1:length(run_scores)
			if run_scores(r) ~= 0
				fprintf(out, '%s\t%d\t%d\t-\t%f\n', chromosome, ...
					pos * granularity, (pos + run_lengths(r)) * granularity, ...
					run_scores(r));
			end
			pos = pos + run_lengths(r);
		end
	end

	fclose(out);
end
