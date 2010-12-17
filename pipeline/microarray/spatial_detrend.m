function samples = spatial_detrend(samples, probes)

probes.XPos = probes.XPos + (1 - min(probes.XPos));
probes.YPos = probes.YPos + (1 - min(probes.YPos));

progress = Progress;

window_size = 7;
window_elems = window_size^2;
radius = (window_size - 1) / 2;

for s = 1:size(samples, 2)
	h = max(probes.YPos) + 2*radius;
	w = max(probes.XPos) + 2*radius;
	f = nan(h, w);
	
	probe_idx = (probes.XPos + radius) * h + probes.YPos + radius;
	
	f(probe_idx) = samples(:, s);
	
	gmed = nanmedian(samples(:, s));

	ignore = find(isnan(f));
	f(ignore) = (mod(ignore, 2) - 0.5) * Inf;
	
	f = medfilt2(f, [window_size window_size]);
	f(ignore) = NaN;
	
	f(f == -Inf) = NaN;
	f(f == Inf) = NaN;
	
	samples(:, s) = samples(:, s) .* (gmed ./ f(probe_idx));
	
	samples(samples(:, s) == Inf, s) = NaN;
	samples(samples(:, s) == -Inf, s) = NaN;
	
	progress.update(s / size(samples, 2));
end

