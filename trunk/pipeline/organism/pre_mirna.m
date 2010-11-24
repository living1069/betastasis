function pre_mirnas = pre_mirna(name)
global organism;
idx = pre_mirna_idx(name);
na = find(isnan(idx));
if length(na) > 0
	error('Could not find pre-miRNA %s.', name{na(1)});
end
pre_mirnas = filter_struct(organism.pre_miRNA, idx);

