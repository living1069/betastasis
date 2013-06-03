function [] = analyze_correlation(x, y)

valid = ~isnan(x) & ~isnan(y);
x = x(valid); y = y(valid);

uniq_x = unique(x);
uniq_y = unique(y);

% Visualization
figure;
if length(uniq_x) <= 8 && length(uniq_y) <= 5
	fprintf('Both groups are categorized, using stacked bar plot...\n');
	stacked = zeros(length(uniq_x), length(uniq_y));
	for i = 1:length(uniq_x)
		for j = 1:length(uniq_y)
			stacked(i, j) = sum(x == uniq_x(i) & y == uniq_y(j));
		end
	end
	for i = 1:length(uniq_x)
		stacked(i, :) = stacked(i, :) / sum(stacked(i, :));
	end
	
	fprintf('Number of samples in groups:');
	for i = 1:length(uniq_x)
		fprintf(' %d', sum(x == uniq_x(i)));
	end
	fprintf('\n');
	
	p = kruskalwallis(y, x);
	fprintf('p = %e (Kruskal-Wallis test)...\n', p);
	
	bar(uniq_x, stacked, 'stacked');
	legend(cellfun(@num2str, num2cell(uniq_y), 'UniformOutput', false));
elseif length(uniq_x) <= 8
	fprintf('X is categorized, using a box plot...\n');
	
	fprintf('Number of samples in groups:');
	for i = 1:length(uniq_x)
		fprintf(' %d', sum(x == uniq_x(i)));
	end
	fprintf('\n');
	
	p = kruskalwallis(y, x);
	fprintf('p = %e (Kruskal-Wallis test)...\n', p);

	boxplot(y, x);
else
	fprintf('Neither is categorized: correlation, scatter plot...\n');
	jx = x;
	if length(x) / length(uniq_x) >= 1.2
		fprintf('Applying jitter to X variable...\n');
		stdev = (max(x) - min(x)) / (length(uniq_x) - 1) / 2;
		jx = x + randn(size(x)) * stdev;
	end

	jy = y;
	if length(y) / length(uniq_y) >= 1.2
		fprintf('Applying jitter to Y variable...\n');
		stdev = (max(y) - min(y)) / (length(uniq_y) - 1) / 2;
		jy = y + randn(size(y)) * stdev;
	end

	scatter(jx, jy, 40, '.');
	
	pearson = corr(x, y, 'type', 'pearson', 'rows', 'pairwise');
	spearman = corr(x, y, 'type', 'spearman', 'rows', 'pairwise');

	fprintf('Pearson correlation: %.2f\n', pearson);
	fprintf('Spearman correlation: %.2f\n', spearman);
end

saveas(gcf, '~/data.pdf');


