function fmatrix = fmatrix_mutations(mutations, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

samples = mutations.Meta.Sample.ID;
if length(samples) ~= length(unique(samples))
	fprintf(1, 'Dataset contains technical replicates. Merging them now...\n');
	test = merge_replicates(test);
end

all_mut_strs = {};
mut_strs = {};
for s = 1:length(mutations.Mutations)
	mut_strs{s} = {};
	for m = 1:length(mutations.Mutations{s}.Gene)
		mut = mutations.Mutations{s};
		mut_strs{s}{m} = sprintf('N:MUT:%s:%s:%d:%d:?', mut.Gene{m}, ...
			chromosomes.Name{mut.Chromosome(m)}, mut.Position(m), ...
			mut.Position(m));
		all_mut_strs{end+1, 1} = mut_strs{s}{m};
	end
end

all_mut_strs = unique(all_mut_strs);

F = length(all_mut_strs);
S = length(mut_strs);

mut_matrix = zeros(F, S);
for s = 1:S
	[~, loc] = ismember(mut_strs{s}, all_mut_strs);
	mut_matrix(loc(loc ~= 0), s) = 1;
end

% Some samples may have no mutation information. Make sure that we place
% only NaN values in those columns.
mut_matrix(:, strcmpi('None', mutations.Method)) = NaN;

fmatrix.Samples = samples;
fmatrix.Features = all_mut_strs;
fmatrix.Data = mut_matrix;
	
fprintf(1, '%d mutation features.\n', F);

