
% CNA_RECURRENCE_TRACK   Calculate a recurrence track of copy number alterations
%
%    CNA_RECURRENCE_TRACK(SEGMENTS, PROBESETS, TRACK_FILE) writes out an IGV
%    track file at path TRACK_FILE, containing information about the recurrence
%    of copy number alterations at different chromosomal loci. The argument
%    PROBESETS must contain CGH probeset information for the aCGH microarray
%    platform in question. The recurrence track is calculated based on segmented
%    copy number data SEGMENTS for each sample pair.
%
%    CNA_RECURRENCE_TRACK(..., 'NormalThreshold', THRESHOLD) specifies the
%    threshold for how high a copy number aberration value must be before
%    being considered as a true aberration. The default value is 0.5, which
%    means that for an autosomal chromosome a copy number of > 2.5 is considered
%    as a true amplification event, and a copy number of < 1.5 as a deletion.

% Author: Matti Annala <matti.annala@tut.fi>

function [] = cna_recurrence_track(segments, cgh_probesets, track_file, ...
	varargin)

global organism;

normal_threshold = 0.5;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'NormalThreshold')
		normal_threshold = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isempty(regexpi(track_file, '.+\.igv'))
	fprintf(1, ['WARNING: A recurrence track should have a .igv suffix for ' ...
	            'optimal visualization in IGV.\n']);
end

cna = cn_seg_expand(segments, cgh_probesets);

S = size(cna, 2);
amp_recurrence = sum(cna > normal_threshold, 2) / S;
del_recurrence = sum(cna < -normal_threshold, 2) / S;

cna_sign = ((amp_recurrence >= del_recurrence) - 0.5) * 2;
recurrence = cna_sign .* max(amp_recurrence, del_recurrence);

fprintf(1, 'Writing CNA recurrence profile as an IGV track...\n');

fid = fopen(track_file, 'W');
fprintf(fid, ...
	'#track maxHeightPixels=500:400:300 graphType=bar viewLimits=--1:1\n');
fprintf(fid, 'Chromosome\tStart\tEnd\tFeature\tCNA_Recurrence\n');

% Find all probes that target the chromosome we're interested in.
for chr = 1:24
	idx = find(cgh_probesets.Chromosome == chr);
	N = length(idx);
	
	offsets = cgh_probesets.Offset(idx);
	chr_recurrence = recurrence(idx);
	
	borders = round(mean([offsets(1:N-1) offsets(2:N)], 2));
	for p = 2:N-1
		fprintf(fid, '%d\t%d\t%d\t-\t%f\n', chr, borders(p-1)+1, borders(p), ...
			chr_recurrence(p));
	end
end

fclose(fid);

