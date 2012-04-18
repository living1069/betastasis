
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

range = [];
show_channels = false;

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Range')
		range = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if strcmpi(varargin{k}, 'ShowChannels')
		show_channels = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

logratios = cgh_to_logratios(samples, refs, probesets, varargin{:});

S = size(logratios, 2);

fid = fopen(track_file, 'W');
fprintf(fid, ...
	'#track maxHeightPixels=500:400:300 graphType=points viewLimits=-2:2\n');
fprintf(fid, 'Chromosome\tStart\tEnd\tFeature');

for s = 1:S, fprintf(fid, '\t%s', samples.meta.sample_id{s}); end
fprintf(fid, '\n');

% Check if the user only wished to write logratio tracks for a small window
% of the genome.
probeset_range = [1 length(probesets.Offset)];
if ~isempty(range)
	tokens = regexpi(range, '^chr(.+?):\s*(\d+)\s*-\s*(\d+)$', 'tokens');
	if length(tokens) ~= 1, error 'Invalid range specified.'; end
	
	token = tokens{1};
	chr = chromosome_sym2num(token{1});
	chr_range = [str2double(token{2}) str2double(token{3})];
	
	if isnan(chr), error 'Invalid chromosome specified.'; end
	
	ps_in_range = find(probesets.Chromosome == chr & ...
		probesets.Offset >= chr_range(1) & probesets.Offset <= chr_range(2));
	probeset_range = [min(ps_in_range) max(ps_in_range)];
end

% Write the logratio tracks into a file.
progress = Progress;
for k = probeset_range(1):probeset_range(2)
	fprintf(fid, '%s\t%d\t%d\t-', ...
		organism.Chromosomes.Name{probesets.Chromosome(k)}, ...
		probesets.Offset(k), probesets.Offset(k));
	fprintf(fid, '\t%f', logratios(k, :));
	fprintf(fid, '\n');
	
	progress.update((k - probeset_range(1)) / ...
		(probeset_range(2) - probeset_range(1)));
end

fclose(fid);

