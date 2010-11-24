
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
	fprintf(1, ['WARNING: Copy number tracks should have a .seg suffix for ' ...
	            'optimal visualization in IGV.\n']);
end

fid = fopen(track_file, 'W');
fprintf(fid, 'Track\tChromosome\tStart\tEnd\tCNA\n');

S = size(segs.Chromosome, 2);

if isfield(segs.Meta, 'Sample') && isfield(segs.Meta.Sample, 'ID')
	track_ids = segs.Meta.Sample.ID;
else
	fprintf(1, ['WARNING: Sample IDs not available. Using surrogate names ' ...
	            'based on sample indices instead.\n']);
	track_ids = {};
	for s = 1:S, track_ids{s} = num2str(s); end
end

% Find all probes that target the chromosome we're interested in.
progress = Progress;
for s = 1:S
	for chr = 1:size(segs.Chromosome, 1)
		chr_segs = segs.Chromosome{chr, s};
		cn = chr_segs.CNA + 2;
		
		for p = 1:length(cn)
			fprintf(fid, '%s\t%d\t%d\t%d\t%f\n', track_ids{s}, chr, ...
				chr_segs.Start(p), chr_segs.End(p), cn(p));
		end
	end
	
	progress.update(s / S);
end

fclose(fid);

