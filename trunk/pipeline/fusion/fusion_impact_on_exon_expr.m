function [] = fusion_impact_on_exon_expr(fusion, rearrangements, exon_expr)

global organism;
genes = organism.Genes;
exons = organism.Exons;

% Discard fusion samples for which we have no exon expression data.
found = ismember(rearrangements.Meta.Sample.ID, exon_expr.Meta.Sample.ID);
if any(~found)
	for s = find(~found)'
		fprintf(1, 'No exon expression data found for fusion sample %s.\n', ...
			rearrangements.Meta.Sample.ID{s});
	end
end

rearrangements = filter_query(rearrangements, found);
fusions = rearrangements.Fusions;
S = length(fusions);

% Check if the given gene expression data matches with our fusion data.
[~, idx] = ismember(rearrangements.Meta.Sample.ID, exon_expr.Meta.Sample.ID);
exon_expr = filter_query(exon_expr, idx);

tokens = regexpi(fusion, '(.+):(.+)', 'tokens');
if length(tokens) ~= 1, error 'Invalid fusion specifier.'; end

token = tokens{1};
gene_5p = token{1}; gene_3p = token{2};

gene_5p = gene_idx(gene_5p);
gene_3p = gene_idx(gene_3p);

if isnan(gene_5p), error 'Unknown 5'' gene.'; end
if isnan(gene_3p), error 'Unknown 3'' gene.'; end

	
% First we build a map of all gene pairs shown to participate in fusions.
gfusion = struct;
gfusion.Genes = [gene_5p, gene_3p];
gfusion.Exons = zeros(0, 2);
gfusion.ReadSamples = {};

for s = 1:S
	for f = 1:length(fusions{s}.ReadCount)
		left_exon = fusions{s}.Exons(f, 1);
		right_exon = fusions{s}.Exons(f, 2);
		
		left_gene = exons.Gene(fusions{s}.Exons(f, 1));
		right_gene = exons.Gene(fusions{s}.Exons(f, 2));
		
		if ~(left_gene == gene_5p && right_gene == gene_3p), continue, end
		
		R = fusions{s}.ReadCount(f);
			
		idx = find(gfusion.Exons(:, 1) == left_exon & ...
			gfusion.Exons(:, 2) == right_exon);
		if isempty(idx)
			gfusion.Exons(end+1, :) = [left_exon, right_exon];
			gfusion.ReadSamples{end+1, 1} = repmat(s, R, 1);
		else
			gfusion.ReadSamples{idx} = [ gfusion.ReadSamples{idx}; ...
				repmat(s, R, 1) ];
		end
	end
end

%if 0
num_reads = cellfun(@length, gfusion.ReadSamples);
if isempty(num_reads)
	fprintf(1, 'Found no samples with both fusion and exon expression data.\n');
	return;
end

[~, largest] = max(num_reads); largest = largest(1);

exon_5p = gfusion.Exons(largest, 1);
exon_3p = gfusion.Exons(largest, 2);

fprintf(1, ['Picked %s[%s]:%s[%s] as the fusion variant with strongest ' ...
	'evidence.\n'], genes.Name{gene_5p}, exons.ID{exon_5p}, ...
	genes.Name{gene_3p}, exons.ID{exon_3p});
	
samples_with_fusion = false(S, 1);
samples_with_fusion(gfusion.ReadSamples{largest}) = true;
%else
%samples_with_fusion = strcmpi(rearrangements.Meta.Sample.ID, 'TCGA-AG-A02X-01');
%find(samples_with_fusion)
%exon_5p = 152497;
%exon_3p = 120801;
%end

exons_5p_gene = find(exons.Gene == gfusion.Genes(1));
exons_3p_gene = find(exons.Gene == gfusion.Genes(2));

[sorted_ids, order] = sort_nat(exons.ID(exons_5p_gene));
exons_5p_5p = exons_5p_gene( ...
	order(1:find(strcmp(exons.ID{exon_5p}, sorted_ids))));
exons_5p_3p = exons_5p_gene( ...
	order(find(strcmp(exons.ID{exon_5p}, sorted_ids))+1:end));
	
[sorted_ids, order] = sort_nat(exons.ID(exons_3p_gene));
exons_3p_5p = exons_3p_gene( ...
	order(1:find(strcmp(exons.ID{exon_3p}, sorted_ids))-1));
exons_3p_3p = exons_3p_gene( ...
	order(find(strcmp(exons.ID{exon_3p}, sorted_ids)):end));

exon_groups = { exons_5p_5p, exons_5p_3p, exons_3p_5p, exons_3p_3p };
exons.ID(exons_5p_5p)
exons.ID(exons_5p_3p)
exons.ID(exons_3p_5p)
exons.ID(exons_3p_3p)



diff_expr = nan(1, 4);
for k = 1:4
	diff_expr(k) = median( ...
		log2(median(exon_expr.Mean(exon_groups{k}, ...
		samples_with_fusion) + 1, 2)) - ...
		log2(median(exon_expr.Mean(exon_groups{k}, ...
		~samples_with_fusion) + 1, 2)));
end

fprintf(1, 'Fusion found in %d / %d samples.\n', ...
	sum(samples_with_fusion), length(samples_with_fusion));



samples_without_fusion = ~samples_with_fusion;
if any(genes.Chromosome([gene_5p gene_3p]) == 24)
	fprintf(1, ['At least one fused gene is in chrY. ' ...
		'Only including males in analysis...\n']);
	
	if ~isfield(exon_expr.Meta, 'Patient') || ...
		~isfield(exon_expr.Meta.Patient, 'Gender')
		return;
	end
	
	fprintf(1, 'Fusion found in %d / %d males.\n', ...
		sum(samples_with_fusion & ...
		strcmp('Male', exon_expr.Meta.Patient.Gender)),...
		sum(strcmp('Male', exon_expr.Meta.Patient.Gender)));
	fprintf(1, 'Fusion found in %d / %d females.\n', ...
		sum(samples_with_fusion & ...
		strcmp('Female', exon_expr.Meta.Patient.Gender)), ...
		sum(strcmp('Female', exon_expr.Meta.Patient.Gender)));
	samples_without_fusion = samples_without_fusion & ...
		strcmp('Male', exon_expr.Meta.Patient.Gender);
end
	
figure('PaperPosition', [.1 .1 .8 .85], 'PaperUnits', 'normalized');
set(gca, 'Units', 'normalized', 'Position', [0.1 0.3 0.8 0.4]);
hold all;

dotsize = 6;
if max(samples_with_fusion, samples_without_fusion) < 20
	dotsize = 10;
end

scatter(1 + .2*sin(2*pi*rand(sum(samples_with_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_5p_5p, samples_with_fusion), 1)), ...
	dotsize, 'filled');
scatter(2 + .2*sin(2*pi*rand(sum(samples_without_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_5p_5p, samples_without_fusion), 1)), ...
	dotsize, 'filled');
scatter(3 + .2*sin(2*pi*rand(sum(samples_with_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_5p_3p, samples_with_fusion), 1)), ...
	dotsize, 'filled');
scatter(4 + .2*sin(2*pi*rand(sum(samples_without_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_5p_3p, samples_without_fusion), 1)), ...
	dotsize, 'filled');
	
scatter(5 + .2*sin(2*pi*rand(sum(samples_with_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_3p_5p, samples_with_fusion), 1)), ...
	dotsize, 'filled');
scatter(6 + .2*sin(2*pi*rand(sum(samples_without_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_3p_5p, samples_without_fusion), 1)), ...
	dotsize, 'filled');
scatter(7 + .2*sin(2*pi*rand(sum(samples_with_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_3p_3p, samples_with_fusion), 1)), ...
	dotsize, 'filled');
scatter(8 + .2*sin(2*pi*rand(sum(samples_without_fusion), 1)), ...
	log2(sum(exon_expr.Mean(exons_3p_3p, samples_without_fusion), 1)), ...
	dotsize, 'filled');
	
xlim([0.5 8.5]); ylabel('Log-2 expression');

set(gca, 'xtick', 1:8, 'xticklabel', { ...
	sprintf('%s-5p (fusion)', genes.Name{gene_5p}), ...
	sprintf('%s-5p (no fusion)', genes.Name{gene_5p}), ...
	sprintf('%s-3p (fusion)', genes.Name{gene_5p}), ...
	sprintf('%s-3p (no fusion)', genes.Name{gene_5p}), ...
	sprintf('%s-5p (fusion)', genes.Name{gene_3p}), ...
	sprintf('%s-5p (no fusion)', genes.Name{gene_3p}), ...
	sprintf('%s-3p (fusion)', genes.Name{gene_3p}), ...
	sprintf('%s-3p (no fusion)', genes.Name{gene_3p})});
	
rotateticklabel(gca, 90);

saveas(gcf, sprintf('%s-%s_scatter.pdf', genes.Name{gene_5p}, ...
	genes.Name{gene_3p}));
	
	
	


