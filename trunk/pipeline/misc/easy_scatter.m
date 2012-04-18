function [] = easy_scatter(a, b, file_prefix, label_y, label_x)

if nargin <= 3, label_y = ''; end
if nargin <= 4, label_x = ''; end

pearson = corr(a, b, 'type', 'pearson', 'rows', 'pairwise')
spearman = corr(a, b, 'type', 'spearman', 'rows', 'pairwise')
	
if nargin > 2
	figure; scatter(a, b, 50, '.');
	xlabel(label_x); ylabel(label_y);
	title(sprintf('Pearson = %.2f         Spearman = %.2f', pearson, spearman));
	saveas(gcf, [file_prefix '.pdf']);
end


