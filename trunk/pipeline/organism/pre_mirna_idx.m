function idx = pre_mirna_idx(name)

global org_;

if ~isfield(org_, 'pre_mirna_name_to_idx')
	global organism;
	org_.pre_mirna_name_to_idx = containers.Map(organism.pre_miRNA.Name, ...
		num2cell(1:length(organism.pre_miRNA.Name)));
end

name_map = org_.pre_mirna_name_to_idx;

if iscellstr(name)
	idx = nan(size(name));
	valid = name_map.isKey(name);
	tmp = cell2mat(name_map.values(name(valid)));
	idx(valid) = tmp;
else
	idx = nan;
	if name_map.isKey(name), idx = name_map(name); end
end

