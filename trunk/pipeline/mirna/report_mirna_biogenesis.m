function [] = report_mirna_biogenesis(test_expr, ref_expr)

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

refinements = zeros(0, 2);
refine_ratios = zeros(0, 2);

for m = 1:length(mirnas.Name)
	[pre, ~] = find(pre_mirnas.Matures == m);
	for p = pre'
		refinements(end+1, :) = [p m];
		refine_ratios(end+1, 1) = median( ...
			test_expr.miRNA(m, :) ./ test_expr.pre_miRNA(p, :));
		refine_ratios(end, 2) = median( ...
			ref_expr.miRNA(m, :) ./ ref_expr.pre_miRNA(p, :));
	end
end

valid = find(all(~isnan(refine_ratios) & refine_ratios ~= Inf & ...
	refine_ratios ~= 0, 2));
[~, order] = sort(abs(log2(refine_ratios(valid, 1)) - ...
	log2(refine_ratios(valid, 2))), 'descend');
order = valid(order);

fprintf(1, 'miRNAs with strongest biogenesis aberration:\n');
for r = order'
	p = refinements(r, 1); m = refinements(r, 2);
	
	if mean(test_expr.miRNA(m, :)) + mean(ref_expr.miRNA(m, :)) < 50 || ...
		mean(test_expr.pre_miRNA(p, :)) + mean(ref_expr.pre_miRNA(p, :)) < 50
		continue;
	end
	
	if 0
	fprintf(1, '%s -> %s: %.2f vs %.2f (test %d -> %d, ref %d -> %d)\n', ...
		pre_mirnas.Name{p}, mirnas.Name{m}, refine_ratios(r, 1), ...
		refine_ratios(r, 2), round(mean(test_expr.pre_miRNA(p, :))), ...
		round(mean(test_expr.miRNA(m, :))), ...
		round(mean(ref_expr.pre_miRNA(p, :))), ...
		round(mean(ref_expr.miRNA(m, :))));
	else
	fprintf(1, '%s\t%s\t%.2f\t%.2f\t%d\t%d\t%d\t%d\n', ...
		pre_mirnas.Name{p}, mirnas.Name{m}, refine_ratios(r, 1), ...
		refine_ratios(r, 2), round(mean(test_expr.pre_miRNA(p, :))), ...
		round(mean(test_expr.miRNA(m, :))), ...
		round(mean(ref_expr.pre_miRNA(p, :))), ...
		round(mean(ref_expr.miRNA(m, :))));
	end
end

