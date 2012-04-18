function reads = reads_for_fusion(fusion_name, fusions)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;
chromosomes = organism.Chromosomes;

fusions = fusions.fusions;

S = length(fusions);

tokens = regexpi(fusion_name, '(.+)[:-](.+)', 'tokens');
if length(tokens) ~= 1, error 'Invalid fusion specifier.'; end

token = tokens{1};
gene_5p_str = token{1}; gene_3p_str = token{2};

gene_5p = gene_idx(gene_5p_str);
gene_3p = gene_idx(gene_3p_str);

if isnan(gene_5p), error('Unknown 5'' gene %s.', gene_5p_str); end
if isnan(gene_3p), error('Unknown 3'' gene %s.', gene_3p_str); end

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

reads = zeros(1, S);
for p = 1:size(gfusion.Exons, 1)
	for s = 1:S
		reads(s) = reads(s) + sum(gfusion.ReadSamples{p} == s);
	end
end

