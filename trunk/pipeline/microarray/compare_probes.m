function probe_changes = compare_probes(probes, ref_probes)
	
seq_to_idx = containers.Map(probes.Sequence, ...
	num2cell(1:length(probes.Sequence)));
ref_seq_to_idx = containers.Map(ref_probes.Sequence, ...
	num2cell(1:length(ref_probes.Sequence)));

N = max(length(probes.Sequence), length(ref_probes.Sequence));
	
added = {};
removed = {};
changed = {};

found = seq_to_idx.isKey(ref_probes.Sequence);
indices = cell2mat(seq_to_idx.values(ref_probes.Sequence(found)));
ref_indices = find(found);

for k = 1:length(indices)
	ref_idx = ref_indices(k);
	idx = indices(k);
	if probes.XPos(idx) ~= ref_probes.XPos(ref_idx) || ...
		probes.YPos(idx) ~= ref_probes.YPos(ref_idx)
		changed{end + 1} = ref_probes.Sequence{ref_idx}; 
	end
end

ref_indices = find(~found);
for k = 1:length(ref_indices)
	ref_idx = ref_indices(k);
	removed{end + 1} = ref_probes.Sequence{ref_idx};
end

indices = find(ref_seq_to_idx.isKey(probes.Sequence));
for k = 1:length(indices)
	idx = indices(k);
	added{end + 1} = probes.Sequence{idx};
end

fprintf(1, 'Probe comparison results:\n');
fprintf(1, '- %d probes added\n', length(added));
fprintf(1, '- %d probes removed\n', length(removed));
fprintf(1, '- %d probes changed position\n', length(changed));

probe_changes = struct( ...
	'Added', { added }, ...
	'Removed', { removed }, ...
	'Changed', { changed });

