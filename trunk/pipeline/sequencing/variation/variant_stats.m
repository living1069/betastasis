
function [] = variant_stats(variants, varargin)

% Calculate a substitution matrix
sub_matrix = zeros(4, 4);
ref = variants.rows.ref_allele; obs = variants.rows.alt_allele;
ref_idx = zeros(size(ref)); obs_idx = zeros(size(obs));

ref_idx(strcmp(ref, 'A')) = 1;
ref_idx(strcmp(ref, 'C')) = 2;
ref_idx(strcmp(ref, 'G')) = 3;
ref_idx(strcmp(ref, 'T')) = 4;
obs_idx(strcmp(obs, 'A')) = 1;
obs_idx(strcmp(obs, 'C')) = 2;
obs_idx(strcmp(obs, 'G')) = 3;
obs_idx(strcmp(obs, 'T')) = 4;

valid = (ref_idx > 0) & (obs_idx > 0);
for k = find(valid)'
	sub_matrix(ref_idx(k), obs_idx(k)) = sub_matrix(ref_idx(k), obs_idx(k)) + 1;
end

sub_matrix
sub_matrix / sum(sum(sub_matrix))

purines = 'AG';
pyrimidines = 'CT';

transversions = { 'A>C', 'C>A', 'A>T', 'T>A', 'G>C', 'C>G', 'G>T', 'T>G' };
transitions = { 'A>G', 'G>A', 'C>T', 'T>C' };

transversions = [ 0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0 ];
transitions =   [ 0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0 ];

num_transversions = sum(sum(sub_matrix .* transversions))
num_transitions = sum(sum(sub_matrix .* transitions))

figure; subplot(221);
pie([num_transversions, num_transitions]);




% Count the number of mutations in each patient.
S = length(variants.meta.sample_id);
total_mutations = nansum(variants.genotype, 1);
total_genotyped = sum(~isnan(variants.genotype), 1);

subplot(222);
title('Total mutations');
pretty_bar(1:S, total_mutations);

subplot(223);
title('Total genotyped');
pretty_bar(1:S, total_genotyped);

subplot(224);
title('Total mutations normalized by total genotyped');
pretty_bar(1:S, total_mutations ./ total_genotyped);

saveas(gcf, '~/mutations_stats.pdf');


