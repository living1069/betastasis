function fmatrix = fmatrix_mirna_expression(expr, varargin)

global organism;
chromosomes = organism.Chromosomes;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

samples = expr.Meta.Sample.ID;
if length(samples) ~= length(unique(samples))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	expr = merge_replicates(expr);
end

valid = find(~any(isnan(expr.Mean), 2));

S = length(expr.Meta.Sample.ID);
F = length(valid);

fmatrix.Samples = expr.Meta.Sample.ID;
fmatrix.Features = cell(F, 1);
fmatrix.Data = log2(expr.Mean(valid, :));

for f = 1:F
	g = valid(f);
	
	[pre, ~] = find(pre_mirnas.Matures == g);
	if isempty(pre) || isnan(pre_mirnas.Chromosome(pre(1)))
		pos_str = '?:?:?:?';
	else
		pre = pre(1);
		pos_str = sprintf('%s:%d:%d:%s', ...
			chromosomes.Name{pre_mirnas.Chromosome(pre)}, ...
			pre_mirnas.Position(pre, 1), pre_mirnas.Position(pre, 2), ...
			pre_mirnas.Strand(pre));
	end

	fmatrix.Features{f} = sprintf('N:EXPR:%s:%s', mirnas.Name{g}, pos_str);
end

fprintf(1, '%d microRNA expression features.\n', F);

