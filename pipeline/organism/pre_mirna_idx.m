function idx = pre_mirna_idx(name)

persistent pre_mirna_name_to_idx;

if isempty(pre_mirna_name_to_idx)
	global organism;
	pre_mirna_name_to_idx = containers.Map(organism.pre_miRNA.Name, ...
		num2cell(1:length(organism.pre_miRNA.Name)));
end

if iscellstr(name)
	idx = nan(size(name));
	valid = pre_mirna_name_to_idx.isKey(name);
	tmp = cell2mat(pre_mirna_name_to_idx.values(name(valid)));
	idx(valid) = tmp;
else
	idx = nan;
	if pre_mirna_name_to_idx.isKey(name)
		idx = pre_mirna_name_to_idx(name);
	end
end

