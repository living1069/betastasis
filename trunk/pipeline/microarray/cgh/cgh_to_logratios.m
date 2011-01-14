
function logratios = cgh_to_logratios(test, ref, probesets, varargin)

global organism;

A = test.Mean;
B = ref.Mean;
S = size(A, 2);

smooth_window_size = 5;
render_histograms = false;
normal_level = nan(S, 1);

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'SmoothWindowSize') || ...
		strcmpi(varargin{k}, 'Smooth')
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

if any(size(A) ~= size(B))
	error 'The sample and reference matrices must have equal dimensions.';
end

if length(normal_level) ~= S
	error 'Length of the zero level argument must equal the amount of samples.';
end

cnv = zeros(length(probesets.ProbeCount), S);
for k = 1:length(probesets.ProbeCount)
	probes = probesets.Probes(k, 1:probesets.ProbeCount(k));
	cnv(k, :) = median(A(probes, :), 1) ./ median(B(probes, :), 1);
end

logratios = log2(cnv);

if smooth_window_size > 0
	for chr = 1:length(organism.Chromosomes.Name)
		idx = find(probesets.Chromosome == chr);
		if isempty(idx), continue, end
			
		a = min(idx); b = max(idx);
		logratios(a:b, :) = medfilt1(logratios(a:b, :), smooth_window_size);
	end
end

for s = 1:S
	if ~isnan(normal_level(s))
		logratios(:, s) = logratios(:, s) - normal_level(s);
		continue;
	end
	
	smoothed = conv2(logratios, ones(9, 1) / 9, 'same');
	
	% Normalize logratios by moving the highest peak to zero on the x-axis.
	bins = -4:0.05:4;
	n = hist(smoothed(:, s), bins);
	
	if render_histograms
		figure; hist(smoothed(:, s), bins);
		xlabel('Probe logratio'); ylabel('Number of probes');
		saveas(gcf, sprintf('%s_hist_%d.pdf', ...
			regexprep(track_file, '\.igv', '', 'ignorecase'), s));
	end

	bins = bins(2:end-1);
	n = n(2:end-1);

	[~, normal_idx] = max(n);
	normal_level(s) = bins(normal_idx);
	
	logratios(:, s) = smoothed(:, s) - normal_level(s);
end

