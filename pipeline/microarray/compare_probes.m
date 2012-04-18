function probe_changes = compare_probes(probes, ref_probes)

% First we check if all probes have maintained their positions.
max_rows = max(max(probes.YPos), max(ref_probes.YPos));
max_cols = max(max(probes.XPos), max(ref_probes.XPos));
seqs = cell(max_rows, max_cols);
ref_seqs = cell(max_rows, max_cols);

seqs(sub2ind(size(seqs), probes.YPos, probes.XPos)) = probes.Sequence;
fprintf(1, 'Test array has %d probes.\n', sum(sum(~cellfun(@isempty, seqs))));

ref_seqs(sub2ind(size(ref_seqs), ref_probes.YPos, ref_probes.XPos)) = ...
	ref_probes.Sequence;
fprintf(1, 'Reference array has %d probes.\n', ...
	sum(sum(~cellfun(@isempty, ref_seqs))));

identical = strcmpi(seqs, ref_seqs);
fprintf(1, '%d identical probes found.\n', sum(sum(identical)));

return;
for y = 1:max_rows
	for x = 1:max_cols
		if identical(y, x) == false
			if isempty(seqs{y,x}) && isempty(ref_seqs{y,x}), continue, end
			fprintf(1, '(%d, %d): %s vs %s\n', x, y, seqs{y,x}, ref_seqs{y,x});
		end
	end
end
