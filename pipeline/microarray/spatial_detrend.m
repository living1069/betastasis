function samples = spatial_detrend(samples, varargin)

window_size = 7;
probes = [];
force_median = [];
log_scale = false;

% For backwards compatibility.
first_option = 1;
if isstruct(varargin{1})
	probes = varargin{1};
	first_option = 2;
end

for k = first_option:2:length(varargin)
	if strcmpi(varargin{k}, 'WindowSize') || strcmpi(varargin{k}, 'Window')
		window_size = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Probes')
		probes = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'ForceMedian')
		force_median = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'LogScale')
		log_scale = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isempty(probes) && isfield(samples, 'Meta')
	probes = platform(samples.Meta.Platform{1}, 'probes');
end

if isfield(samples, 'Mean'), samples = samples.Mean; end
	
if ~log_scale
	if any(samples(:, 1) < 0)
		error(['Found negative values. Use the option "LogScale" if your ' ...
			'data is in logarithmic scale.']);
	end
end

probes.XPos = probes.XPos + (1 - min(probes.XPos));
probes.YPos = probes.YPos + (1 - min(probes.YPos));

progress = Progress;

window_elems = window_size^2;
radius = (window_size - 1) / 2;

for s = 1:size(samples, 2)
	h = max(probes.YPos) + 2*radius;
	w = max(probes.XPos) + 2*radius;
	f = nan(h, w);
	
	probe_idx = (probes.XPos + radius) * h + probes.YPos + radius;
	
	f(probe_idx) = samples(:, s);
	
	if isempty(force_median)
		gmed = nanmedian(samples(:, s));
	else
		gmed = force_median;
	end

	ignore = find(isnan(f));
	f(ignore) = (mod(ignore, 2) - 0.5) * Inf;
	
	f = medfilt2(f, [window_size window_size]);
	f(ignore) = NaN;
	
	f(f == -Inf) = NaN;
	f(f == Inf) = NaN;
	
	if log_scale
		samples(:, s) = samples(:, s) + (gmed - f(probe_idx));
	else
		samples(:, s) = samples(:, s) .* (gmed ./ f(probe_idx));
	end
	
	samples(samples(:, s) == Inf, s) = NaN;
	samples(samples(:, s) == -Inf, s) = NaN;
	
	progress.update(s / size(samples, 2));
end

