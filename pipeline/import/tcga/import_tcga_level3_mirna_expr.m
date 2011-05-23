function meta = import_tcga_level3_mirna_expr(filename, dataset)

global organism;
mirnas = organism.miRNA;

fid = fopen(filename);
data = textscan(fid, '%s %s %f', 'Headerlines', 1, ...
	'Delimiter', '\t', 'ReturnOnError', 0, 'TreatAsEmpty', 'NA');
fclose(fid);

samples = data{1};
mirna = data{2};
mirna_expr = data{3};

expr.Meta.Type = 'miRNA expression';
expr.Meta.Sample.ID = unique(samples);
S = length(expr.Meta.Sample.ID);
expr.Meta.Platform = repmat({'Agilent Human miRNA 8x15K'}, S, 1);

expr.Mean = nan(length(mirnas.Name), S);

sample_to_idx = containers.Map(expr.Meta.Sample.ID, ...
	num2cell(1:S));

for s = 1:S
	lines = find(strcmp(expr.Meta.Sample.ID{s}, samples));
	
	midx = mirna_idx(mirna(lines));
	valid = ~isnan(midx);
	
	expr.Mean(midx(valid), s) = mirna_expr(lines(valid));
end

expr.Mean = 2.^expr.Mean;

fprintf(1, 'Expression values not found for %d miRNA.\n', ...
	sum(any(isnan(expr.Mean), 2)));
	
create_dataset(dataset, expr);
augment_meta_tcga(dataset);




