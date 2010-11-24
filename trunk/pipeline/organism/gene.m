function genes = gene(name)
global organism;
idx = gene_idx(name);
na = find(isnan(idx));
if length(na) > 0
	error('Could not find gene %s.', name{na(1)});
end

genes = filter_struct(organism.Genes, idx);

