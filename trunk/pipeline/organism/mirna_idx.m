function idx = mirna_idx(name)

persistent mirna_name_to_idx;

if isempty(mirna_name_to_idx)
	global organism;
	mirna_name_to_idx = containers.Map(organism.miRNA.Name, ...
		num2cell(1:length(organism.miRNA.Name)));
end

if iscellstr(name)
	idx = nan(size(name));
	valid = mirna_name_to_idx.isKey(name);
	tmp = cell2mat(mirna_name_to_idx.values(name(valid)));
	idx(valid) = tmp;
else
	idx = nan;
	if mirna_name_to_idx.isKey(name)
		idx = mirna_name_to_idx(name);
	end
end

