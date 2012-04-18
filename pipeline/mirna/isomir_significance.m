function [] = isomir_significance(isomirs_test, isomirs_ref, report_dir)

global organism;
mirnas = organism.miRNA;

abs_threshold = 20;

sum_isomirs_test = zeros(length(mirnas.Name), length(isomirs_test{1, 1}.Count));
sum_isomirs_ref = zeros(length(mirnas.Name), length(isomirs_ref{1, 1}.Count));

for m = 1:size(isomirs_test, 1)
	for s = 1:size(isomirs_test, 2)
		isomirs_test{m, s}.Count(isnan(isomirs_test{m, s}.Count)) = 0;
		sum_isomirs_test(m, :) = sum_isomirs_test(m, :) + ...
			isomirs_test{m, s}.Count';
	end
end
for m = 1:size(isomirs_ref, 1)
	for s = 1:size(isomirs_ref, 2)
		isomirs_ref{m, s}.Count(isnan(isomirs_ref{m, s}.Count)) = 0;
		sum_isomirs_ref(m, :) = sum_isomirs_ref(m, :) + ...
			isomirs_ref{m, s}.Count';
	end
end

test_mirna_avg_len = zeros(length(mirnas.Name), 1);
ref_mirna_avg_len = zeros(length(mirnas.Name), 1);

total_test_mirna_reads = sum(sum_isomirs_test, 2);
total_ref_mirna_reads = sum(sum_isomirs_ref, 2);

for m = 1:length(mirnas.Name)
	test_mirna_avg_len(m) = ...
		nansum(isomirs_test{m, 1}.End .* sum_isomirs_test(m, :)) / ...
		total_test_mirna_reads(m);
	ref_mirna_avg_len(m) = ...
		nansum(isomirs_ref{m, 1}.End .* sum_isomirs_ref(m, :)) / ...
		total_ref_mirna_reads(m);
end

%[test_mirna_avg_len, ref_mirna_avg_len]

diff_avg_len = test_mirna_avg_len - ref_mirna_avg_len;
[~, order] = sort(abs(diff_avg_len), 'descend');

[~, ~] = mkdir(report_dir);

fprintf(1, 'Strongest differential miRNA size polymorphisms:\n');
shown = 0;
for k = 1:length(mirnas.Name)
	if shown == 20, break, end
		
	idx = order(k);
	if isnan(diff_avg_len(idx)), continue, end
	if any(mirnas.Name{idx} == '*'), continue, end
		
	if total_test_mirna_reads(idx) < abs_threshold || ...
		total_ref_mirna_reads(idx) < abs_threshold
		continue;
	end
	
	fprintf(1, '- %s: %f (avg. counts %d, %d)\n', mirnas.Name{idx}, ...
		diff_avg_len(idx), ...
		round(total_test_mirna_reads(idx) / size(isomirs_test, 2)), ...
		round(total_ref_mirna_reads(idx) / size(isomirs_ref, 2)));
	shown = shown + 1;
	
	safe_name = strrep(mirnas.Name{idx}, '*', 'star');
	isomir_hist(isomirs_test, isomirs_ref, mirnas.Name{idx}, ...
		[report_dir '/' safe_name '.pdf']);
end

