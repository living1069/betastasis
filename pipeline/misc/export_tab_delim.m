function [] = export_tab_delim(expr, filename)

global organism;

drop_nan = true;

if strcmpi(expr.Meta.Type, 'Gene expression')
	out = fopen(filename, 'W');
	progress = Progress;
	
	fprintf(out, 'NAME');
	for k = 1:length(expr.Meta.Sample.ID)
		fprintf(out, '\t%s', expr.Meta.Sample.ID{k});
	end
	fprintf(out, '\n');
	
	if drop_nan
		valid_genes = find(~any(isnan(expr.Mean), 2));
	else
		valid_genes = 1:length(organism.Genes.Name);
	end
	
	for g = valid_genes'
		fprintf(out, '%s', organism.Genes.Name{g});
		fprintf(out, '\t%f', expr.Mean(g, :));
		fprintf(out, '\n');
		progress.update(g / length(organism.Genes.Name));
	end
else
	error(['Exporting to tab delimited format is not supported for this ' ...
	       'type of data.']);
end

fclose(out);

