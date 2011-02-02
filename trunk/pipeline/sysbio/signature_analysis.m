
% Author: Matti Annala <matti.annala@tut.fi>

function [] = signature_analysis(test, ref, links, varargin)

global organism;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Significance')
		significance = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

test_log_expr = log2(test.Mean);
ref_log_expr = log2(ref.Mean);

regulators = unique(links.Genes(:, 1));

progress = Progress;

min_S = min(size(test_log_expr, 2), size(ref_log_expr, 2))
%max_dim = min_S - 1;
%max_dim = round(min_S / 2);

max_dims = [4 9 14 19 24 29];
dist_score = nan(length(organism.Genes.Name), length(max_dims));

for d = 1:length(max_dims)
	max_dim = max_dims(d);

	for r = 1:length(regulators)
		g = regulators(r);
		
		% We have to restrict the subnet to the min_S genes with the highest
		% MI score. Otherwise the sample covariance matrix becomes
		% singular (non-invertible).
		sublinks = find(links.Genes(:, 1) == g);
		[~, order] = sort(links.MI(sublinks), 'descend');
		
		subnet = links.Genes(sublinks(order), 2);
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
		
		progress.update(((d-1) + r / length(regulators)) / length(max_dims));
	end
end

dist_score = median(dist_score, 2);

valid = find(~isnan(dist_score));

[~, order] = sort(dist_score(valid), 'descend');
order = valid(order);

fprintf(1, 'Genes with highest signature delta:\n');
for k = 1:30
	fprintf(1, '- %s (distance %e)\n', organism.Genes.Name{order(k)}, ...
		dist_score(order(k)));
end

