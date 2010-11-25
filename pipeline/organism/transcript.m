function transcripts = transcript(name)
global organism;
if isnumeric(name) && isvector(name)
	idx = name;
	for k = 1:length(idx)
		if idx < 0 || idx > length(organism.Transcripts.Name)
			error('Transcript with index %d does not exist.', idx(k));
		end
	end
else
	if ischar(name)
		name = { name };
	end
	idx = transcript_idx(name);
	na = find(isnan(idx));
	if length(na) > 0
		error('Could not find transcript %s.', name{na(1)});
	end
end

transcripts = filter_struct(organism.Transcripts, idx);

