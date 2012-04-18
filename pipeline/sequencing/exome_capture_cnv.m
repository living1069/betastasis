function [] = exome_capture_cnv(bam_test, bam_ref, range, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

read_len = 76;
resolution = 1000;
bandwidth = 2 * resolution;
informative_density = 20;
plot_height = 0.2;

ylimits = [-3 3];
cbs_alpha = 0.01;
show_genes = {};
show_exons = '';

for k = 1:2:length(varargin)
	if regexpi(varargin{k}, 'YLim')
		ylimits = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if strcmpi(varargin{k}, 'ShowGenes')
		show_genes = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if strcmpi(varargin{k}, 'ShowExons')
		show_exons = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end

tokens = regexpi(range, '^(chr.+):\s*(\d+)\s*-\s*(\d+)$', 'tokens');
if length(tokens) ~= 1, error 'Invalid range.'; end

token = tokens{1};
chr = chromosome_sym2num(token{1});
range = str2double({token{2}, token{3}});

test_al = bamread(bam_test, token{1}, range);
ref_al = bamread(bam_ref, token{1}, range);

test_pos = nan(length(test_al), 1);
ref_pos = nan(length(ref_al), 1);
for k = 1:length(test_pos), test_pos(k) = test_al(k).Position; end
for k = 1:length(ref_pos), ref_pos(k) = ref_al(k).Position; end

%fid = fopen(bam_test);
%test_data = textscan(fid, '%*s %d %s %d %*s %*[^\n]', 'Delimiter', '\t', ...
%	'ReturnOnError', false);
%fclose(fid);

range = [range(1)-read_len, range(2)+read_len]

x = range(1):resolution:range(2);

test_density = ksdensity(test_pos + read_len/2, x, ...
	'kernel', 'box', 'width', bandwidth);
ref_density = ksdensity(ref_pos + read_len/2, x, ...
	'kernel', 'box', 'width', bandwidth);

% This sucks!
test_density = test_density / min(test_density(test_density > 0));
ref_density = ref_density / min(ref_density(ref_density > 0));

valid = (test_density > informative_density) | ...
	(ref_density > informative_density);
	
test_density = test_density + 1;
ref_density = ref_density + 1;

logratio = log2(test_density ./ ref_density);

if 0

figure('PaperUnits', 'normalized', 'PaperPosition', [0.0 0.75 1.0 plot_height]);
%subplot('Position', [0.0 0.75 1.0 plot_height]);
hold all;

if ~isempty(show_genes)
	for k = 1:length(show_genes)
		g = show_genes{k};
		if ischar(g)
			gene_pos = genes.Position(gene_idx(g));
		elseif isnumeric(g)
			gene_pos = g;
		end
		fill(gene_pos([1 1 2 2]), ylimits([1 2 2 1]), [.9 .9 .9]);
	end
end

if ~isempty(show_exons)
	ex = find(organism.Exons.Gene == gene_idx(show_exons));
	exon_ids = organism.Exons.ID(ex);
	exon_pos = organism.Exons.Position(ex, :);
	for e = 1:size(exon_pos, 1)
		fill(exon_pos(e, [1 1 2 2]), ylimits([1 2 2 1]), [.9 .9 .9]);
		text(mean(exon_pos(e, :)), ...
			ylimits(2) + mod(e, 2) * (ylimits(2) - ylimits(1)) / 20, ...
			exon_ids{e}, 'HorizontalAlignment', 'center');
	end
end

scatter(x(valid), logratio(valid), 5, [0 0 0], 'filled');

cbs_in = [repmat(chr, sum(valid), 1), x(valid)', logratio(valid)'];
seg = cghcbs(cbs_in, 'StoppingRule', true, 'Alpha', cbs_alpha);
seg = seg.SegmentData
line([seg.Start'; seg.End'], repmat(seg.Mean', 2, 1), 'Color', 'r', ...
	'LineWidth', 2);
	
xlabel(sprintf('chr%s offset', chromosomes.Name{chr}));
ylabel('Log-2 ratio');
if ~isempty(ylimits), ylim(ylimits); end
xlim(range);
saveas(gcf, 'exome_capture_cnv.pdf');
end



figure('PaperUnits', 'normalized', 'PaperPosition', [0.0 0.75 1.0 plot_height]);
scatter(x, test_density, 5, [0 0 0], 'filled'); xlim(range);
saveas(gcf, 'exome_capture_test_density.pdf');


figure('PaperUnits', 'normalized', 'PaperPosition', [0.0 0.75 1.0 plot_height]);
scatter(x, ref_density, 5, [0 0 0], 'filled'); xlim(range);
saveas(gcf, 'exome_capture_ref_density.pdf');



