function [] = fusion_impact_on_expr(fusion, rearrangements, expr)

global organism;
genes = organism.Genes;
exons = organism.Exons;

fusions = rearrangements.Fusions;
S = length(fusions);

% Check if the given gene expression data matches with our fusion data.
[found, idx] = ismember(rearrangements.Meta.Sample.ID, expr.Meta.Sample.ID);
if any(~found)
	error 'Expression data does not match with rearrangement samples.';
end
expr = filter_query(expr, idx);

tokens = regexpi(fusion, '(.+):(.+)', 'tokens');
if length(tokens) ~= 1, error 'Invalid fusion specifier.'; end

token = tokens{1};
gene_5p = token{1}; gene_3p = token{2};

gene_5p = gene_idx(gene_5p);
gene_3p = gene_idx(gene_3p);

if isnan(gene_5p), error 'Unknown 5'' gene.'; end
if isnan(gene_3p), error 'Unknown 3'' gene.'; end

samples_with_fusion = false(S, 1);
	
for s = 1:S
	for f = 1:length(fusions{s}.ReadCount)
		left_gene = exons.Gene(fusions{s}.Exons(f, 1));
		right_gene = exons.Gene(fusions{s}.Exons(f, 2));
		
		if gene_5p == left_gene && gene_3p == right_gene
			samples_with_fusion(s) = true;
		end
	end
end

fprintf(1, 'Fusion found in %d / %d samples.\n', sum(samples_with_fusion), ...
	length(samples_with_fusion));
fprintf(1, 'Fusion found in %d / %d males.\n', ...
	sum(samples_with_fusion & strcmp('Male', expr.Meta.Patient.Gender)), ...
	sum(strcmp('Male', expr.Meta.Patient.Gender)));
fprintf(1, 'Fusion found in %d / %d females.\n', ...
	sum(samples_with_fusion & strcmp('Female', expr.Meta.Patient.Gender)), ...
	sum(strcmp('Female', expr.Meta.Patient.Gender)));

box_expr = log2([expr.Mean([gene_5p gene_5p], :)'+1, ...
	expr.Mean([gene_3p gene_3p], :)'+1]);
box_expr(samples_with_fusion, [1 3]) = NaN;
box_expr(~samples_with_fusion, [2 4]) = NaN;

figure; boxplot(box_expr, 'labels', { ...
		sprintf('%s (fusion)', genes.Name{gene_5p}), ...
		sprintf('%s (no fusion)', genes.Name{gene_5p}), ...
		sprintf('%s (fusion)', genes.Name{gene_3p}), ...
		sprintf('%s (no fusion)', genes.Name{gene_3p}) }, ...
	'labelorientation', 'inline');

saveas(gcf, sprintf('~/%s-%s_boxplot.pdf', genes.Name{gene_5p}, ...
	genes.Name{gene_3p}));

%samples_without_fusion = ~samples_with_fusion & strcmp('Male', expr.Meta.Patient.Gender);
samples_without_fusion = ~samples_with_fusion;
	
figure; hold all;
scatter(1 + .3*sin(2*pi*rand(sum(samples_with_fusion), 1)), ...
	log2(expr.Mean(gene_5p, samples_with_fusion)));
scatter(2 + .3*sin(2*pi*rand(sum(samples_without_fusion), 1)), ...
	log2(expr.Mean(gene_5p, samples_without_fusion)));
scatter(3 + .3*sin(2*pi*rand(sum(samples_with_fusion), 1)), ...
	log2(expr.Mean(gene_3p, samples_with_fusion)));
scatter(4 + .3*sin(2*pi*rand(sum(samples_without_fusion), 1)), ...
	log2(expr.Mean(gene_3p, samples_without_fusion)));

set(gca, 'xtick', 1:4, 'xticklabel', { ...
	sprintf('%s (fusion)', genes.Name{gene_5p}), ...
	sprintf('%s (no fusion)', genes.Name{gene_5p}), ...
	sprintf('%s (fusion)', genes.Name{gene_3p}), ...
	sprintf('%s (no fusion)', genes.Name{gene_3p}) });

saveas(gcf, sprintf('~/%s-%s_scatter.pdf', genes.Name{gene_5p}, ...
	genes.Name{gene_3p}));

