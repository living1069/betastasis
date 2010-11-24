function probes = intersect_probes(probes_to_filter, ref_probes)

useq = unique(ref_probes.Sequence);
ref_seq_map = containers.Map(useq, num2cell(zeros(length(useq), 1)));

N = length(probes_to_filter.Sequence);

probes_to_keep = false(N, 1);

fprintf(1, 'Filtering probes...\n');
fprintf(1, 'Progress: 00%%');

progress = 0;
	
for k = 1:N
	if ref_seq_map.isKey(probes_to_filter.Sequence{k})
		probes_to_keep(k) = 1;
	end
	
	if floor(k / N * 100) > progress
		progress = floor(k / N * 100);
		fprintf(1, '\b\b\b%02d%%', progress);
	end
end

fprintf(1, '\n');

probes = struct('XPos', probes_to_filter.XPos(probes_to_keep), ...
				'YPos', probes_to_filter.YPos(probes_to_keep), ...
				'Sequence', { probes_to_filter.Sequence(probes_to_keep) });

fprintf(1, 'Kept %d probes out of %d.\n', sum(probes_to_keep), N);
fprintf(1, 'Matching was done against %d reference probes.\n', ...
	length(ref_probes.Sequence));
