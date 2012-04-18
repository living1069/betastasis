
% Author: Matti Annala <matti.annala@tut.fi>

function survival_genes = survival_associated_genes(expr, varargin)

global organism;

genes_to_report = 100;

log_expr = log2(expr.mean);
[G, S] = size(log_expr);

feature_chr = nan(G, 1);

if rx(expr.meta.type, 'gene expression')
	features = expr.rows.gene_symbol;
	gidx = gene_idx(features);
	valid = ~isnan(gidx);
	feature_chr(valid) = organism.Genes.Chromosome(gidx(valid));
elseif rx(expr.meta.type, 'mirna expression')
	features = expr.rows.mirna_symbol;
	gidx = mirna_idx(features);
	valid = ~isnan(gidx);
	[~, pos] = ismember(gidx(valid), organism.pre_miRNA.Matures(:));
	[pre_mirna_idx, ~] = ind2sub(size(organism.pre_miRNA.Matures), pos);
	feature_chr(valid) = organism.pre_miRNA.Chromosome(pre_mirna_idx);
else
	error 'Expression dataset is of unrecognized type.';
end

check = 'hsa-miR-21$';
fprintf('At position %d.\n', find(rx(features, check)));

valid = valid & all(~isnan(log_expr), 2);
fprintf('%d / %d (%.1f%%) features with expression for all samples.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);
fprintf('Present after validation: %d\n', valid(rx(features, check)));

valid = valid & mad(log_expr, 1, 2) > 0.2;
fprintf('%d / %d (%.1f%%) features after discarding low variance features.\n',...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);
fprintf('Present after validation: %d\n', valid(rx(features, check)));

valid = valid & feature_chr <= 22;
fprintf('%d / %d (%.1f%%) features after discarding sex chromosomes.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);
fprintf('Present after validation: %d\n', valid(rx(features, check)));
	
valid = find(valid);


clusters = cell(G, 1);
vrs = nan(G, 1);
logrank_p = nan(G, 1);







fprintf('Calculating bimodality scores...\n');

kmeans_start = [min(log_expr, [], 2), max(log_expr, [], 2)];

for g = valid'
	clusters{g} = kmeans(log_expr(g, :)', 2, 'start', kmeans_start(g, :)');
	
	tss = sum((log_expr(g, :) - mean(log_expr(g, :))).^2);
	wss = sum((log_expr(g, clusters{g} == 1) - ...
		mean(log_expr(g, clusters{g} == 1))).^2) + ...
		sum((log_expr(g, clusters{g} == 2) - ...
		mean(log_expr(g, clusters{g} == 2))).^2);
	vrs(g) = wss / tss;
end

% Run a logrank test on the survival profiles for the top bimodal genes.
[~, order_vrs] = sort(vrs);

for g = order_vrs(1:min(genes_to_report, length(valid)))'
	low = clusters{g} == 1;
	high = clusters{g} == 2;
	
	if mean(log_expr(g, low)) > mean(log_expr(g, high))
		error 'WTF centroids flipped?';
	end
	
	threshold = mean([max(log_expr(g, low)), min(log_expr(g, high))]);
	
	figure; subplot(211); hold all;
	[f, xi] = ksdensity(log_expr(g, :));
	plot(xi(xi < threshold), f(xi < threshold), ...
		'Color', 'blue', 'LineWidth', 2);
	plot(xi(xi >= threshold), f(xi >= threshold), ...
		'Color', 'red', 'LineWidth', 2);
	
	xlabel(sprintf('Log-2 expression of %s', features{g}));
	ylabel('Probability density estimate');
	
	subplot(212);
	logrank_p(g) = logrank(...
		expr.meta.survival_time(clusters{g} == 1), ...
		expr.meta.survival_time_censored(clusters{g} == 1), ...
		expr.meta.survival_time(clusters{g} == 2), ...
		expr.meta.survival_time_censored(clusters{g} == 2));
		
	legend('Low expression', 'High expression', 'Censored');
	
	filename = sprintf('%s.pdf', features{g});
	filename = strrep(filename, '*', 'star');
	saveas(gcf, filename);
end

% Then we rank the top genes according to their effect on survival.
[~, order_logrank] = sort(logrank_p);

fid = fopen('top_bimodal_survival_features.txt', 'W');
fprintf(fid, 'GENE\tLOGRANK_TEST_P_VAL_LOG10\tBIMODALITY_VRS\n');

for g = order_logrank(1:min(genes_to_report, length(valid)))'
	fprintf(fid, '%s\t%.2f\t%.3f\n', features{g}, -log10(logrank_p(g)), ...
		vrs(g));
end

fclose(fid);

