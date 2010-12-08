
% CGH_LOGRATIO_TRACK   Visualize CGH logratios as an IGV track
%
%    CGH_LOGRATIO_TRACK(TEST, REF, PROBESETS, TRACK_FILE) constructs an IGV
%    track TRACK_FILE that contains CGH logratios between the paired samples
%    in TEST and REF. CGH probesets must be provided as the PROBESETS argument.
%
%    CGH_LOGRATIO_TRACK(..., 'Smooth', WS) tells the function to
%    first median filter the logratios using a window of size WS. This reduces
%    the noise level in the logratio track. Set to 0 for no smoothing.
%    Default is to smooth using a window of size 5.
%
%    CGH_LOGRATIO_TRACK(..., 'RenderHistograms', true) tells the function to
%    render a probe logratio histogram for each sample. The histograms are
%    stored as PDFs with the same filename prefix as TRACK_FILE.
%
%    CGH_LOGRATIO_TRACK(..., 'NormalLevel', NLEVEL) tells the function to not
%    use automatic logratio level normalization, and instead uses the normal
%    levels given by the user in NLEVEL. NLEVEL must be a vector of S elements,
%    where S is the number of samples in TEST and REF. A value of NaN in NLEVEL
%    tells the function to use automatic level normalization, while any other
%    value tells the function to use the given logratio level for normal
%    (unaberrated) regions.
%
%    See also PAIRED_SAMPLES, CN_SEG_TO_TRACK.

% Author: Matti Annala <matti.annala@tut.fi>

function [] = cgh_logratio_track(samples, refs, probesets, track_file, varargin)

global organism;

logratios = cgh_to_logratios(samples, refs, probesets, varargin{:});

S = size(logratios, 2);

fid = fopen(track_file, 'W');
fprintf(fid, ...
	'#track maxHeightPixels=500:400:300 graphType=points viewLimits=-2:2\n');
fprintf(fid, 'Chromosome\tStart\tEnd\tFeature');

if isfield(samples.Meta, 'Sample') && isfield(samples.Meta.Sample, 'ID')
	for s = 1:S
		fprintf(fid, '\t%s', samples.Meta.Sample.ID{s});
	end
else
	fprintf(1, ['WARNING: Sample IDs not available. Using surrogate names ' ...
	            'based on sample indices instead.\n']);
	for s = 1:S
		fprintf(fid, '\t%d', s);
	end
end
fprintf(fid, '\n');

progress = Progress;
for k = 1:size(logratios, 1)
	fprintf(fid, '%s\t%d\t%d\t-', ...
		organism.Chromosomes.Name{probesets.Chromosome(k)}, ...
		probesets.Offset(k), probesets.Offset(k));
	fprintf(fid, '\t%f', logratios(k, :));
	fprintf(fid, '\n');
	
	progress.update(k / size(logratios, 1));
end

fclose(fid);

