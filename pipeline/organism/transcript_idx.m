function idx = transcript_idx(name)

persistent transcript_name_to_idx;

if isempty(transcript_name_to_idx)
	global organism;
	transcript_name_to_idx = containers.Map(organism.Transcripts.Name, ...
		num2cell(1:length(organism.Transcripts.Name)));
end

if iscellstr(name)
	idx = nan(size(name));
	valid = transcript_name_to_idx.isKey(name);
	tmp = cell2mat(transcript_name_to_idx.values(name(valid)));
	idx(valid) = tmp;
else
	idx = nan;
	if transcript_name_to_idx.isKey(name)
		idx = transcript_name_to_idx(name);
	end
end

