function mirnas = mirna(name)
global organism;
idx = mirna_idx(name);
na = find(isnan(idx));
if length(na) > 0
	error('Could not find miRNA %s.', name{na(1)});
end

mirnas = filter_struct(organism.miRNA, idx);

