function [] = isomir_hist(isomirs_test, isomirs_ref, mirna, report_file)

global organism;
mirnas = organism.miRNA;

m = mirna_idx(mirna);

test_isomir_dist = zeros(1, length(isomirs_test{1,1}.Count));
ref_isomir_dist = zeros(1, length(isomirs_ref{1,1}.Count));

for s = 1:size(isomirs_test, 2)
	isomirs_test{m, s}.Count(isnan(isomirs_test{m, s}.Count)) = 0;
	test_isomir_dist(1, :) = test_isomir_dist(1, :) + isomirs_test{m, s}.Count';
end
for s = 1:size(isomirs_ref, 2)
	isomirs_ref{m, s}.Count(isnan(isomirs_ref{m, s}.Count)) = 0;
	ref_isomir_dist(1, :) = ref_isomir_dist(1, :) + isomirs_ref{m, s}.Count';
end

figure('PaperOrientation', 'landscape', 'PaperPosition', [.1 .1 .8 .8], ...
	'PaperUnits', 'normalized');

subplot(121);
bar(isomirs_test{m, 1}.End, test_isomir_dist / sum(test_isomir_dist));
xlabel('Isomir length'); ylabel('Ratio of reads');
title([mirna ' isomirs in CRC'], 'Interpreter', 'none');
	
subplot(122);
bar(isomirs_ref{m, 1}.End, ref_isomir_dist / sum(ref_isomir_dist));
xlabel('Isomir length'); ylabel('Ratio of reads');
title([mirna ' isomirs in endometrial'], 'Interpreter', 'none');

saveas(gcf, report_file);

