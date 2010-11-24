function idx = gene_idx(name)

persistent gene_name_to_idx;

if isempty(gene_name_to_idx)
	global organism;
	gene_name_to_idx = containers.Map(organism.Genes.Name, ...
		num2cell(1:length(organism.Genes.Name)));
end

if iscellstr(name)
	idx = nan(size(name));
	valid = gene_name_to_idx.isKey(name);
	tmp = cell2mat(gene_name_to_idx.values(name(valid)));
	idx(valid) = tmp;
else
	idx = nan;
	if gene_name_to_idx.isKey(name)
		idx = gene_name_to_idx(name);
	end
end

