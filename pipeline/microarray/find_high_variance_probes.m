function bad_probes = find_high_variance_probes(samples, probesets)

% Calculate linear coefficient for the relationship between probe intensity
% and probe standard deviation.

probe_means = mean(samples.Mean, 2);
probe_stdevs = std(samples.Mean, 0, 2);

correlation = corr([probe_means, probe_stdevs]);
fprintf(1, ['Pearson correlation between standard deviation and probe ' ...
            'intensity is %f.\n'], correlation(2, 1));

theta = [ones(length(probe_means), 1), probe_means] \ probe_stdevs;

normalized_probe_stdevs = probe_stdevs - (theta(1) + theta(2) * probe_means);
%normalized_probe_stdevs = normalized_probe_stdevs - ...
%	mean(normalized_probe_stdevs);

threshold = 14 * std(normalized_probe_stdevs);

bad_mask = normalized_probe_stdevs > threshold;

figure; hold all;
scatter(probe_means(~bad_mask), normalized_probe_stdevs(~bad_mask));
scatter(probe_means(bad_mask), normalized_probe_stdevs(bad_mask));
saveas(gcf, '~/mean_stdev_fit.png');

fprintf(1, 'Found %d bad probes.\n', sum(bad_mask));

bad_probes = find(bad_mask);
%figure;
%N = length(normalized_probe_stdevs(~bad_mask));
%scatter(sort(normalized_probe_stdevs(~bad_mask) / ...
%	std(normalized_probe_stdevs(~bad_mask))), norminv((1:N) / (N + 1))');
%saveas(gcf, '~/qqplot.png');

