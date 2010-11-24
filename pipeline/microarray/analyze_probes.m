function [] = analyze_probes(probes)

unique_sequences = unique(probes.Sequence);

seq_to_num = containers.Map(unique_sequences, ...
	num2cell(zeros(1, length(unique_sequences))));

for p = 1:length(probes.Sequence)
	seq_to_num(probes.Sequence{p}) = seq_to_num(probes.Sequence{p}) + 1;
end

counts = cell2mat(seq_to_num.values(unique_sequences));

fprintf(1, 'Unique probes: %d (%.1f%%)\n', sum(counts == 1), ...
	sum(counts == 1) / length(unique_sequences));

most_redundant = find(counts == max(counts));
fprintf(1, 'Most redundant probes with %d copies:\n', max(counts));
for k = 1:length(most_redundant)
	fprintf(1, '- %s\n', unique_sequences{most_redundant(k)});
end

figure;
hist(counts, 0:max(counts)); xlim([0 max(counts)]);
saveas(gcf, '~/probe_hist.png');
