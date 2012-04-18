
function logratios = cgh_to_logratios(test, ref, probesets, varargin)

global organism;

A = test.mean;
if isempty(ref), ref.mean = ones(size(A)); end
B = ref.mean;

S = size(A, 2);

smooth_window_size = 5;
render_histograms = false;
normal_level = nan(S, 1);
%detrend = false;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'SmoothWindowSize') || ...
		strcmpi(varargin{k}, 'Smooth')
		smooth_window_size = varargin{k+1}; continue;
	end
	
	if strcmpi(varargin{k}, 'Detrend')
		detrend = varargin{k+1}; continue;
	end
	
	if strcmpi(varargin{k}, 'RenderHistograms')
		render_histograms = varargin{k+1}; continue;
	end
	
	if strcmpi(varargin{k}, 'NormalLevel')
		normal_level = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if any(size(A) ~= size(B))
	error 'The sample and reference matrices must have equal dimensions.';
end

if length(normal_level) ~= S
	error 'Length of the zero level argument must equal the amount of samples.';
end

logratios = log2(A ./ B);

%if detrend
%	probes = platform(test.Meta.Platform{1}, 'probes');
%	logratios = spatial_detrend(logratios, 'Probes', probes, ...
%		'ForceMedian', 0, 'WindowSize', 9, 'LogScale', true);
%end

% Transform the probe logratios into probeset logratios. If an array has probes
% with identical sequences, average the logratios across them.
if all(probesets.ProbeCount == 1)
	ps_logr = logratios(probesets.Probes(:, 1), :);
else
	ps_logr = zeros(length(probesets.ProbeCount), S);
	for k = 1:length(probesets.ProbeCount)
		probes = probesets.Probes(k, 1:probesets.ProbeCount(k));
		ps_logr(k, :) = median(logratios(probes, :), 1);
	end
end

logratios = ps_logr;

% Replace negative infinities with the smallest non-infinite value.
logratios(logratios == -Inf) = NaN;
logratios(isnan(logratios)) = nanmin(nanmin(logratios));

if smooth_window_size > 0
	for chr = 1:length(organism.Chromosomes.Name)
		idx = find(probesets.Chromosome == chr);
		if isempty(idx), continue, end
			
		a = min(idx); b = max(idx);
		logratios(a:b, :) = medfilt1(logratios(a:b, :), smooth_window_size);
	end
end

%if detrend, return, end


	
smoothed = conv2(logratios, ones(15, 1) / 15, 'same');

for s = 1:S
	if ~isnan(normal_level(s))
		logratios(:, s) = logratios(:, s) - normal_level(s);
		continue;
	end
	
	% Normalize logratios by moving the highest peak to zero on the x-axis.
	bins = quantile(smoothed(:, s), 0.2)-0.1: ...
		0.02:quantile(smoothed(:, s)+0.1, 0.8);
	n = hist(smoothed(:, s), bins);
	
	if render_histograms
		figure; hist(smoothed(:, s), bins);
		xlabel('Probe logratio'); ylabel('Number of probes');
		saveas(gcf, sprintf('hist_%d.pdf', s));
	end

	bins = bins(2:end-1);
	n = n(2:end-1);
	
	[~, normal_idx] = max(n);
	normal_level(s) = bins(normal_idx);
	
	logratios(:, s) = logratios(:, s) - normal_level(s);
end

