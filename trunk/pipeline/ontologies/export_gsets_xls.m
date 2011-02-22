function [] = export_gsets_xls(gsets, xls_file)

global organism;
genes = organism.Genes;

fields = fieldnames(gsets);
use_function = true;
if ~isempty(find(strcmp('Genes', fields)))
	use_function = false;
end

fid = fopen(xls_file, 'W');
for k = 1:length(gsets.Name)
	if isempty(gsets.Name{k}), continue, end
	
	fprintf(fid, '%s', gsets.Name{k});
	
	if use_function == false
		gset = gsets.Genes{k};
	else
		gset = gsets.genes(k);
	end
	
	for g = 1:length(gset)
		fprintf(fid, '\t%s', genes.Name{gset(g)});
	end
	fprintf(fid, '\n');
end
fclose(fid);

