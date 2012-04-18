function fusions = fusion_convert_old(fusions)

S = length(fusions.Fusions);

for s = 1:S
	% Convert old format read sequences to the new format.
	if size(fusions.Fusions{s}.ReadSequences, 2) > 1
		new_seq = cell(size(fusions.Fusions{s}.ReadSequences, 1), 1);
		for k = 1:size(fusions.Fusions{s}.ReadSequences, 1)
			new_seq{k} = fusions.Fusions{s}.ReadSequences(k, ...
				1:fusions.Fusions{s}.ReadCount(k))';
		end
		fusions.Fusions{s}.ReadSequences = new_seq;
		
		new_joffsets = cell(size(fusions.Fusions{s}.ReadSequences, 1), 1);
		for k = 1:size(fusions.Fusions{s}.ReadSequences, 1)
			new_joffsets{k} = fusions.Fusions{s}.JunctionOffsets(k, ...
				1:fusions.Fusions{s}.ReadCount(k))';
		end
		fusions.Fusions{s}.JunctionOffsets = new_joffsets;
	end
end


