
% Author: Matti Annala <matti.annala@tut.fi>

function dist_score = signature_analysis(test, ref, links, varargin)

global organism;

metric = 'bhatta';
report = '';

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Metric')
		metric = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Report')
		report = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

test_log_expr = log2(test.Mean);
ref_log_expr = log2(ref.Mean);

regulators = unique(links.Genes(:, 1));

if ~isempty(report)
	[~, ~] = mkdir(report);
	[~, ~] = mkdir([report '/mdscale']);
	[~, ~] = mkdir([report '/target_expr']);
end


if regexpi(metric, 'bhatta')

	progress = Progress;
	
	min_S = min(size(test_log_expr, 2), size(ref_log_expr, 2))

	max_dims = [4 9 14 19 24 29];
	%max_dims = 29;
	dist_score = nan(length(organism.Genes.Name), length(max_dims));

	for d = 1:length(max_dims)
		max_dim = max_dims(d);

		for r = 1:length(regulators)
			g = regulators(r);
			
			% We have to restrict the subnet to the min_S genes with the highest
			% MI score. Otherwise the sample covariance matrix becomes
			% singular (non-invertible).
			if 1
				sublinks = find(links.Genes(:, 1) == g);
				[~, order] = sort(links.MI(sublinks), 'descend');
				subnet = links.Genes(sublinks(order), 2);
			else
				subnet = links.Genes(links.Genes(:, 1) == g, 2);
				[~, ttest_p] = ttest2(test_log_expr(subnet, :)', ...
					ref_log_expr(subnet, :)');
				[~, order] = sort(ttest_p, 'ascend');
				subnet = subnet(order);
			end
			
			if length(subnet) > max_dim
				subnet = subnet(1:max_dim);
			end
			subnet = [subnet; g];
			
			sub_test = test_log_expr(subnet, :);
			sub_ref = ref_log_expr(subnet, :);
			
			mu_test = mean(sub_test, 2);
			mu_ref = mean(sub_ref, 2);
			
			sigma_test = cov(sub_test');
			sigma_ref = cov(sub_ref');
			
			sigma = (sigma_test + sigma_ref) / 2;
			
			dist_score(g, d) = 1/8 * (mu_test - mu_ref)' * inv(sigma) * ...
				(mu_test - mu_ref) + ...
				1/2 * log(det(sigma) / sqrt(det(sigma_test) * det(sigma_ref)));
			
			progress.update(((d-1) + r / length(regulators)) / ...
				length(max_dims));
		end
	end

	dist_score = median(dist_score, 2);

elseif regexpi(metric, 'zscore')
	
	dist_score = nan(length(organism.Genes.Name), 1);
	
	for r = 1:length(regulators)
		g = regulators(r);
		
		subnet = links.Genes(links.Genes(:, 1) == g, 2);
		
		sub_test = nanmean(test_log_expr(subnet, :), 2);
		sub_ref = nanmean(ref_log_expr(subnet, :), 2);
		
		dist_score(g) = nanmean(sub_test - sub_ref);
	end
	
else
	error 'Specified distance metric is not supported.';
end
	

valid = find(~isnan(dist_score));
[~, order] = sort(abs(dist_score(valid)), 'descend');
order = valid(order);

order(50) = gene_idx('AR');

fprintf(1, 'Genes with highest signature delta:\n');
for k = 1:50
	fprintf(1, '- %s (distance %.2f)\n', organism.Genes.Name{order(k)}, ...
		dist_score(order(k)));
end

if ~isempty(report)
	figure; hist(dist_score, 100);
	xlabel('Bhattacharyya distance'); ylabel('Number of TFs');
	saveas(gcf, [report '/dist_score.pdf']);
	
	max_dim = 29;
	score_order = order;
	
	for k = 1:50
		g = score_order(k);
		sublinks = find(links.Genes(:, 1) == g);
		[~, order] = sort(links.MI(sublinks), 'descend');
		subnet = links.Genes(sublinks(order), 2);

		if length(subnet) > max_dim
			subnet = subnet(1:max_dim);
		end
		subnet = [g; subnet];

		% Render figures showing how well the samples are grouped according
		% to the gene signature.
		test_idx = 1:size(test_log_expr, 2);
		ref_idx = (1:size(ref_log_expr, 2)) + size(test_log_expr, 2);
		
		sub_expr = [test_log_expr(subnet, :)'; ref_log_expr(subnet, :)'];
		D = squareform(pdist(sub_expr));

		points = mdscale(D, 2);
		figure; hold all;
		scatter(points(test_idx, 1), points(test_idx, 2));
		scatter(points(ref_idx, 1), points(ref_idx, 2));
		saveas(gcf, [report '/mdscale/' organism.Genes.Name{g} '.pdf']);
		
		% Render figures showing how the expression levels of the correlated
		% genes differ between the two sample groups.
		subnet = subnet(end:-1:1);
		figure; barh(mean(test_log_expr(subnet, :), 2) - ...
			mean(ref_log_expr(subnet, :), 2));
		set(gca, 'ytick', 1:length(subnet), ...
			'yticklabel', organism.Genes.Name(subnet));
		set(gca, 'FontSize', 8);
		title(sprintf('Gene neighborhood signature for %s', ...
			organism.Genes.Name{g}));
		xlabel('Differential log2 expression');
		saveas(gcf, [report '/target_expr/' organism.Genes.Name{g} '.pdf']);
	end
end


