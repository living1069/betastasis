
% CGH_LOGRATIO_TRACK   Visualize CGH logratios as an IGV track
%
%    CGH_LOGRATIO_TRACK(TEST, REF, PROBESETS, TRACK_FILE) constructs an IGV
%    track TRACK_FILE that contains CGH logratios between the paired samples
%    in TEST and REF. CGH probesets must be provided as the PROBESETS argument.
%
%    CGH_LOGRATIO_TRACK(..., 'SmoothWindowSize', WS) tells the function to
%    first median filter the logratios using a window of size WS. This reduces
%    the noise level in the logratio track. Set to 0 for no smoothing.
%    Default is to smooth using a window of size 7.
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

function [] = cgh_logratio_track(samples, refs, probesets, track_file, varargin)

global organism;

A = samples.Mean;
B = refs.Mean;
S = size(A, 2);

smooth_window_size = 7;
render_histograms = false;
normal_level = nan(S, 1);

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'SmoothWindowSize')
		smooth_window_size = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'RenderHistograms')
		render_histograms = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'NormalLevel')
		normal_level = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end


if length(normal_level) ~= S
	error 'Length of the zero level argument must equal the amount of samples.';
end

if isempty(regexpi(track_file, '.+\.igv'))
	fprintf(1, 'WARNING: CGH logratio tracks should have a .igv suffix.\n');
end

cnv = zeros(length(probesets.ProbeCount), S);
for k = 1:length(probesets.ProbeCount)
	probes = probesets.Probes(k, 1:probesets.ProbeCount(k));
	cnv(k, :) = median(A(probes, :), 1) ./ median(B(probes, :), 1);
end

logratios = log2(cnv);

if smooth_window_size > 0
	for chr = 1:24
		idx = find(probesets.Chromosome == chr);
		a = min(idx); b = max(idx);
		
		logratios(a:b, :) = medfilt2(logratios(a:b, :), [smooth_window_size 1]);
	end
end

for s = 1:S
	if ~isnan(normal_level(s))
		logratios(:, s) = logratios(:, s) - normal_level(s);
		continue;
	end
	
	% Normalize logratios by moving the highest peak to zero on the x-axis.
	bins = -4:0.05:4;
	n = hist(logratios(:, s), bins);
	
	if render_histograms
		figure; hist(logratios(:, s), bins);
		xlabel('Probe logratio'); ylabel('Number of probes');
		saveas(gcf, sprintf('%s_hist_%d.pdf', ...
			regexprep(track_file, '\.igv', '', 'ignorecase'), s));
	end

	bins = bins(2:end-1);
	n = n(2:end-1);

	[~, normal_idx] = max(n);
	normal_level(s) = bins(normal_idx);
	
	logratios(:, s) = logratios(:, s) - normal_level(s);
end

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



if 0
	fprintf(1, 'Calculating Lomb periodogram...\n');
	idx = find(probesets.Chromosome == 1);
	a = min(idx); b = max(idx);
	fastlomb(logratios(a:b, 1), probesets.Offset(a:b), gcf);
	saveas(gcf, '~/lomb.pdf');
end

