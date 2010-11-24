function gender = cgh_detect_gender(samples, probesets)

N = length(probesets.ProbeCount);
S = size(samples, 2);

chr_intensity = zeros(24, S);
for chr = 1:24
	chr_probesets = (probesets.Chromosome == chr);
	probes = probesets.Probes(chr_probesets, :);
	probes = probes(:);
	probes = probes(probes ~= 0);
	chr_intensity(chr, :) = median(samples(probes, :), 1);
end

chr_x_ratio = chr_intensity(23, :) ./ median(chr_intensity(1:22, :), 1);
chr_y_ratio = chr_intensity(24, :) ./ median(chr_intensity(1:22, :), 1);

% Finally we have a formula for this.
%manliness = min(chr_y_ratio * 2, 1);
manliness = 0.5 * min(chr_y_ratio * 2, 1) + 0.5 * min(2 * max(1 - chr_x_ratio, 0), 1);

gender = repmat({'Male'}, S, 1);
for k = 1:S
	if manliness(k) < 0.5
		gender{k} = 'Female';
	end
end

