
% Author: Matti Annala <matti.annala@tut.fi>

function [] = cgh_plot_logratios(test, refs, probesets, range, ...
	varargin)

global organism;
genes = organism.Genes;

show_genes = {};
show_exons = '';
segments = [];
file_prefix = 'cgh';
plot_height = 0.2;
user_ylimits = [];
cbs_alpha = 0.01;

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if rx(varargin{k}, 'prefix')
		file_prefix = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if rx(varargin{k}, 'segment')
		segments = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if rx(varargin{k}, 'ylim')
		user_ylimits = varargin{k+1};
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
varargin = varargin(~drop_args);

logratios = cgh_to_logratios(test, refs, probesets, 'Smooth', 0, ...
	varargin{:});
	
S = size(logratios, 2);

% Parse the genome region for which the user wishes to show the logratios.
tokens = regexpi(range, '^chr(.+?):\s*(\d+)\s*-\s*(\d+)$', 'tokens');
if length(tokens) ~= 1, error 'Invalid range specified.'; end
	
token = tokens{1};
chr = chromosome_sym2num(token{1});
view_range = [str2double(token{2}) str2double(token{3})];

if isnan(chr), error 'Invalid chromosome specified.'; end

	
	


ps_in_view = find(probesets.Chromosome == chr & ...
	probesets.Offset >= view_range(1) & probesets.Offset <= view_range(2));
	
for s = 1:S
	figure('PaperUnits', 'normalized', ...
		'PaperPosition', [0.0 0.75 1.0 plot_height]);
	hold all;
	
	if isempty(user_ylimits)
		ylimits = quantile(logratios(ps_in_view, s), [.01 .99]);
		dy = (ylimits(2) - ylimits(1)) / 5;
		ylimits = [ylimits(1)-dy, ylimits(2)+dy]
	else
		ylimits = user_ylimits;
	end

	
	if ~isempty(show_genes)
		for g = gene_idx(show_genes)
			fill(genes.Position(g, [1 1 2 2]), ylimits([1 2 2 1]), [.9 .9 .9]);
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

	scatter(probesets.Offset(ps_in_view), logratios(ps_in_view, s), ...
		5, [0 0 0], 'filled');
	
	if ~isempty(segments)
		seg = segments.chromosome{chr, s}
		line([seg.Start'; seg.End'], repmat(seg.Mean', 2, 1), 'Color', 'r', ...
			'LineWidth', 2);
	end
		
	xlabel(sprintf('chr%s offset', organism.Chromosomes.Name{chr}));
	ylabel('Log-2 ratio');
	ylim(ylimits);
	xlim([view_range(1) view_range(2)]);
	saveas(gcf, [file_prefix '_' test.meta.sample_id{s} '.pdf']);
end


