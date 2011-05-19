function fmatrix = fmatrix_diff_expression(test, ref, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = test.Meta.Sample.ID;
if length(samples) ~= length(unique(samples))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	test = merge_replicates(test);
end

test_expr = log2(test.Mean);
ref_expr = log2(ref.Mean);

valid = find(~any(isnan(ref_expr), 2));

S = length(test.Meta.Sample.ID);
F = length(valid);

fmatrix.Samples = test.Meta.Sample.ID;
fmatrix.Features = cell(F, 1);
fmatrix.Data = nan(F, S);

for f = 1:F
	g = valid(f);
	
	data = test_expr(g, :) - median(ref_expr(g, :));
	
	if isnan(genes.Chromosome(g))
		pos_str = '?:?:?:?';
	else
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{genes.Chromosome(g)}, ...
			genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
	end
	
	fmatrix.Features{f} = sprintf('N:DIFF_EXPR:%s:%s', genes.Name{g}, pos_str);
	fmatrix.Data(f, :) = data;
end

fprintf(1, '%d differential expression features.\n', F);

