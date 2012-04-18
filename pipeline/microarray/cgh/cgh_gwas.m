
% Author: Matti Annala <matti.annala@tut.fi>

function p = cgh_gwas(raw, confounders, response, probesets, varargin)

global organism;

if length(response) ~= size(raw.Mean, 2)
	error 'Size of response variable does not match with number of samples.';
end

if any(isnan(response))
	fprintf(1, ['Response variable contains %d / %d missing observations, ' ...
		'ignoring them.\n'], sum(isnan(response)), length(response));
	
	valid = ~isnan(response);
	confounders = confounders(valid, :);
	response = response(valid);
	raw = filter_query(raw, valid);
end

S = size(raw.Mean, 2);

iqr_plot(raw, '~/iqr_plot_before_normalization.pdf');


% Perform global median normalization.
medians = nanmedian(raw.Mean, 1);
med_medians = median(medians);
for s = 1:S
	raw.Mean(:, s) = raw.Mean(:, s) / (medians(s) / med_medians);
end


% Drop outliers.
qtiles = quantile(raw.Mean, [.05 .25 .5 .75 .95]);
[~, outliers] = deleteoutliers(qtiles(5, :), 0.05);

fprintf(1, '%d / %d (%.1f%%) samples were declared as outliers.\n', ...
	length(outliers), S, length(outliers) / S * 100);

iqr_plot(raw, '~/iqr_plot_after_normalization.pdf', 'Order', 'IQR', ...
	'Mark', outliers);

confounders = confounders(setdiff(1:S, outliers), :);
response = response(setdiff(1:S, outliers));
raw = filter_query(raw, setdiff(1:S, outliers));
S = size(raw.Mean, 2);

iqr_plot(raw, '~/iqr_plot_after_outlier_removal.pdf', 'Order', 'IQR');





% If a single genomic locus has more than one probe assigned to it, we combine
% their signals by taking a median across the redundant probe intensities.
ps_raw = summarize_probesets(raw, probesets);


% Use only the top quartile of probes that exhibit the highest variation.
kept_probes = size(ps_raw, 1);
%kept_probes = round(size(ps_raw, 1) / 4);

score = mad(ps_raw, 1, 2) ./ median(ps_raw, 2);
[~, order] = sort(score, 'descend');
ps_raw(order(kept_probes+1:end), :) = NaN;


% Divide the data into two groups based on the (Boolean) response variable.
A = ps_raw(:, response == 1);
B = ps_raw(:, response == 2);



%figure; ksdensity(median(ps_raw, 2)); saveas(gcf, '~/foo_range.pdf');

%figure; ksdensity(ps_raw(900, :)); saveas(gcf, '~/foo_900.pdf');
%figure; ksdensity(ps_raw(950, :)); saveas(gcf, '~/foo_950.pdf');
%figure; ksdensity(ps_raw(1000, :)); saveas(gcf, '~/foo_1000.pdf');









[~, p] = ttest2(A', B', 0.05, 'both', 'unequal');
[~, order] = sort(p, 'ascend');

fdr = 0.05;

num_significant = 0;
for k = kept_probes:-1:1
	if p(order(k)) < k * fdr / length(order)
		num_significant = k;
	end
end

%bonferroni_sig = 0.05 / size(A, 1);
%num_significant = sum(p <= bonferroni_sig);

fprintf(1, '%d probes significant after Benjamini-Hochberg FDR control.\n', ...
	num_significant);

order = order(1:1000);

% For the probes that we found significant based on t-test, we now perform
% logistic regression to determine their odds ratio.
lr_stats = nan(size(ps_raw, 1), 2);

fprintf(1, 'TOP GWAS REGIONS (single probe):\n');
fprintf(1, 'Chromosome\tPosition\tP (t-test)\tP (LR)\tOdds ratio\n');
for k = order
	[coeff, ~, stats] = mnrfit([ps_raw(k, :)', confounders], response);
	lr_stats(k, :) = [coeff(1), stats.p(1)];
end

[~, order] = sort(lr_stats(:, 2)', 'ascend');

for k = order(1:25)
	control = nanmedian(B(k, :));
	
	% Render a scatter plot of the probe intensity in the two response
	% categories.
	figure; hold all;
	[f, xi] = ksdensity(A(k, :) / control);
	plot(xi, f, 'r', 'LineWidth', 2);
	[f, xi] = ksdensity(B(k, :) / control);
	plot(xi, f, 'k', 'LineWidth', 2);
	
	scatter(A(k, :) / control, 5 + randn(1, size(A, 2)) / 5, 10, 'r', 'filled');
	scatter(B(k, :) / control, 5 + randn(1, size(B, 2)) / 5, 10, 'k', 'filled');
	
	saveas(gcf, sprintf('~/chr%s:%d.pdf', ...
		organism.Chromosomes.Name{probesets.Chromosome(k)}, ...
		probesets.Offset(k)));
	
	fprintf(1, 'chr%s\t%d\t%.3e\t%.3e\t%.2f\n', ...
		organism.Chromosomes.Name{probesets.Chromosome(k)}, ...
		probesets.Offset(k), p(k), lr_stats(k, 2), exp(lr_stats(k, 1)));
end






% MEDIAN 3
p_3 = 10 .^ -medfilt1(-log10(p), 3);
[~, order] = sort(p_3, 'ascend');

fprintf(1, 'TOP GWAS REGIONS (median 3):\n');
for k = order(1:50)
	fprintf(1, 'chr%s - %d: p-value %.3e\n', ...
		organism.Chromosomes.Name{probesets.Chromosome(k)}, ...
		probesets.Offset(k), p(k));
end

% MEDIAN 5
p_5 = 10 .^ -medfilt1(-log10(p), 5);
[~, order] = sort(p_5, 'ascend');

fprintf(1, 'TOP GWAS REGIONS (median 5):\n');
for k = order(1:50)
	fprintf(1, 'chr%s - %d: p-value %.3e\n', ...
		organism.Chromosomes.Name{probesets.Chromosome(k)}, ...
		probesets.Offset(k), p(k));
end





track_file = '~/foo.igv';
fid = fopen(track_file, 'W');
fprintf(fid, ...
	'#track maxHeightPixels=500:400:300 graphType=points viewLimits=0:15\n');
fprintf(fid, 'Chromosome\tStart\tEnd\tFeature\tP-value\n');

for k = 1:length(probesets.Chromosome)
	fprintf(fid, '%s\t%d\t%d\t-\t%f\n', ...
		organism.Chromosomes.Name{probesets.Chromosome(k)}, ...
		probesets.Offset(k), probesets.Offset(k), -log10(p(k)));
end

fclose(fid);





function raw = normalize_zero_level(raw)

for s = 1:size(raw, 2)
	% Normalize logratios by moving the highest peak to zero on the x-axis.
	qtile = quantile(raw(:, s), [0.2 0.8]);
	bins = qtile(1):0.05:qtile(2);
	n = hist(raw(:, s), bins);
	
	bins = bins(2:end-1);
	n = n(2:end-1);
	
	[~, normal_idx] = max(n);
	normal_level(s) = bins(normal_idx);
	
	raw(:, s) = raw(:, s) - normal_level(s);
end





function ps_raw = summarize_probesets(raw, probesets)

if all(probesets.ProbeCount == 1)
	ps_raw = raw.Mean(probesets.Probes(:, 1), :);
else
	ps_raw = zeros(length(probesets.ProbeCount), size(raw.Mean, 2));
	for k = 1:length(probesets.ProbeCount)
		probes = probesets.Probes(k, 1:probesets.ProbeCount(k));
		ps_raw(k, :) = median(raw.Mean(probes, :), 1);
	end
end

