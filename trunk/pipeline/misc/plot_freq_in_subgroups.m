function [] = plot_freq_in_subgroups(groups, group_names, ...
	pos_label, pos_samples, neg_label, neg_samples)

features = nan(length(groups), 2);
norm_features = nan(length(groups), 2);
for f = 1:length(groups)
	features(f, 1) = sum(groups{f} & pos_samples);
	features(f, 2) = sum(groups{f} & neg_samples);
	norm_features(f, :) = features(f, :) / sum(features(f, :));
end

figure('PaperOrientation', 'landscape', ...
	'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
bar(100 * norm_features, 'stacked');
legend(sprintf('%s (n = %d)', pos_label, sum(pos_samples)), ...
	sprintf('%s (n = %d)', neg_label, sum(neg_samples)));
ylim([0 120]); ylabel('Percentage of patients');
set(gca, 'XTick', 1:size(features, 1), 'XTickLabel', group_names);
rotateticklabel(gca, 45);

% Add subgroup sample counts.
for f = 1:length(groups)
	ypos = 103;
	if length(groups) > 10 && mod(f, 2) == 0
		ypos = 106;
	end
	text(f, ypos, sprintf('%d / %d', features(f, 1), sum(features(f, :))), ...
		'HorizontalAlignment', 'center');
end

saveas(gcf, '~/barchart.pdf');

