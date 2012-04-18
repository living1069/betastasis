function fmatrix = fmatrix_fusions(fusions, varargin)

global organism;
genes = organism.Genes;

sample_ids = exon_expr.Meta.Sample.ID;
if length(sample_ids) ~= length(unique(sample_ids))
	error 'Dataset contains technical replicates. Cannot continue...';
end

joint_fusions = filter_fusions(fusions, varargin);

S = length(fusions.Fusions);
F = length(joint_fusions);

fusion_matrix = zeros(F, S);
for f = 1:F
	for p = 1:size(gfusion.Exons, 1)
		fusion_matrix(f, gfusion.ReadSamples{p}) = 1;
	end
end

fmatrix.Samples = fusions.Meta.Sample.ID;
fmatrix.Features = cell(F, 1);
fmatrix.Data = fusion_matrix;

for f = 1:F
	fmatrix.Features{f} = sprintf('N:FUSION:%s-%s', ...
		genes.Name{joint_fusions{f}.Genes(1)}, ...
		genes.Name{joint_fusions{f}.Genes(2)});
end

fprintf(1, '%d fusion gene features.\n', F);

