
% CN_SEG_TO_TRACK     Export segmented copy number samples as IGV tracks
%
%    CN_SEG_TO_TRACK(SEGS, TRACK_FILE) exports the segmented copy number samples
%    SEGS as individual tracks in TRACK_FILE.
%
%    It is recommended that the filename in TRACK_FILE end in a .seg suffix,
%    since this is the standard suffix used in IGV for segmented data files.

% Author: Matti Annala <matti.annala@tut.fi>

function [] = cn_seg_to_track(segs, track_file)

if isempty(regexpi(track_file, '.+\.seg'))
	fprintf(['WARNING: Copy number tracks should have a .seg suffix for ' ...
		'optimal visualization in IGV.\n']);
end

fid = fopen(track_file, 'W');
fprintf(fid, 'Track\tChromosome\tStart\tEnd\tCNA\n');

S = size(segs.chromosome, 2);

track_ids = segs.meta.sample_id;

% Find all probes that target the chromosome we're interested in.
progress = Progress;
for s = 1:S
	for chr = 1:size(segs.chromosome, 1)
		chr_segs = segs.chromosome{chr, s};
		lr = chr_segs.logratio;
		
		for p = 1:length(lr)
			fprintf(fid, '%s\t%d\t%d\t%d\t%f\n', track_ids{s}, chr, ...
				chr_segs.start(p), chr_segs.end(p), lr(p));
		end
	end
	
	progress.update(s / S);
end

fclose(fid);

