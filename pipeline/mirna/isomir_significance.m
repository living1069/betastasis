function [] = isomir_significance(isomirs_test, isomirs_ref)

global organism;
mirnas = organism.miRNA;

test_mirna_avg_len = zeros(size(isomirs_test));
ref_mirna_avg_len = zeros(size(isomirs_ref));

test_total_count = zeros(length(mirnas.Name), 1);
ref_total_count = zeros(length(mirnas.Name), 1);

for m = 1:size(test_mirna_avg_len, 1)
	for s = 1:size(test_mirna_avg_len, 2)
		total = nansum(isomirs_test{m, s}.Count);
		test_mirna_avg_len(m, s) = nansum(isomirs_test{m, s}.End' .* ...
			isomirs_test{m, s}.Count / total);
		test_total_count(m) = test_total_count(m) + total;
	end
end

for m = 1:size(ref_mirna_avg_len, 1)
	for s = 1:size(ref_mirna_avg_len, 2)
		total = nansum(isomirs_ref{m, s}.Count);
		ref_mirna_avg_len(m, s) = nansum(isomirs_ref{m, s}.End' .* ...
			isomirs_ref{m, s}.Count / total);
		ref_total_count(m) = ref_total_count(m) + total;
	end
end

test_mirna_avg_len(test_mirna_avg_len == 0) = NaN;
ref_mirna_avg_len(ref_mirna_avg_len == 0) = NaN;

diff_avg_len = nanmean(test_mirna_avg_len, 2) - nanmean(ref_mirna_avg_len, 2)
[~, order] = sort(abs(diff_avg_len), 'descend');

fprintf(1, 'Strongest differential miRNA size polymorphisms:\n');
shown = 0;
for k = 1:length(mirnas.Name)
	if shown == 50, break, end
	idx = order(k);
	
	if test_total_count(idx) < 20 || ref_total_count(idx) < 20
		continue;
	end
	
	fprintf(1, '- %s: %f (counts %d, %d)\n', mirnas.Name{idx}, ...
		diff_avg_len(idx), test_total_count(idx), ref_total_count(idx));
	shown = shown + 1;
end

