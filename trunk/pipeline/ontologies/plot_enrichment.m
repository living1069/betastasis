
% Author: Matti Annala <matti.annala@tut.fi>

function pval = plot_enrichment(values, positive, range)

gsea_ylim = [-.3 .3];

if isnumeric(positive)
	positive = positive(~isnan(positive));
	tmp = positive;
	positive = false(size(values));
	positive(tmp) = true;
end

% Allow inputs to be column or row vectors.
values = values(:);
positive = positive(:);

valid = find(~isnan(values));
values = values(valid);
positive = positive(valid);

permutations = 100000;

[values, L] = sort(values, 'descend');
hit = positive(L);
Nh = sum(hit);

Phit = cumsum(hit) / Nh;
Pmis = cumsum(~hit) / (length(L) - Nh);

[~, idx] = max(abs(Phit - Pmis));
escore = Phit(idx) - Pmis(idx)

figure; subplot('position', [0.05 0.56 0.9 0.3]); hold all;
%figure; subplot('position', [0.35 0.56 0.3 0.3]); hold all;
plot(Phit - Pmis);
line([0 length(hit)], [0 0], 'Color', [0 0 0]);
ylim(gsea_ylim);
set(gca, 'XTick', []);

if nargin < 3
	range = [min(values) 0 max(values)];
end

norm = values;
norm(values >= range(2)) = (values(values >= range(2)) - range(2)) ./ ...
	range(3) - range(2);
norm(values < range(2)) = -(values(values < range(2)) - range(2)) ./ ...
	range(1) - range(2);
norm(norm > 1) = 1;
norm(norm < -1) = -1;

subplot('position', [0.05 0.5 0.9 0.05]);
grad = zeros(1, length(Phit), 3);
for k = 1:length(norm)
	if hit(k)
		grad(1, k, :) = [0 0 0];
	elseif norm(k) >= 0
		grad(1, k, :) = [1 1-norm(k) 1-norm(k)];
	else
		grad(1, k, :) = [1+norm(k) 1+norm(k) 1];
	end
end

image(grad);
lines = line([find(hit)'; find(hit)'], repmat([0; 2], 1, sum(hit)));
set(lines, 'Color', [0 0 0]);
set(gca, 'XTick', [], 'YTick', []);
xlim([0 length(hit)]);

subplot('position', [0.05 0 0.9 0.48]);

bars = area(abs(values)); xlim([0 length(values)]);
set(bars, 'LineStyle', 'none', 'FaceColor', [.8 .8 .8]);
set(gca, 'XTick', []);

saveas(gcf, '~/gsea.pdf');






perm_escores = zeros(1, permutations);

Nh = sum(hit);
N = length(hit);

for p = 1:permutations
	hit = false(N, 1);
	r = randperm(N);
	hit(r(1:Nh)) = true;
	
	Phit = cumsum(hit) / Nh;
	Pmis = cumsum(~hit) / (length(L) - Nh);

	[~, idx] = max(abs(Phit - Pmis));
	perm_escores(p) = Phit(idx) - Pmis(idx);
end

if escore >= 0
	pval = sum(perm_escores > escore) / length(perm_escores);
else
	pval = sum(perm_escores < escore) / length(perm_escores);
end
pval = max(pval, 1 / permutations)




%fprintf(1, 'P-value = %.4f\t\tE-score = %.2f\t\t%s = %f\n', ...
%	enriched.Pval(idx), enriched.Escore(idx), ratio_label, ...
%	nanmean(diff(enriched.Genes{idx})));


