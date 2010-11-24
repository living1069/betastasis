function [] = export_tab_delim(expr, filename)

global organism;

if strcmpi(expr.Meta.Type, 'Gene expression')
	out = fopen(filename, 'W');
	progress = Progress;
	
	fprintf(out, 'NAME');
	for k = 1:length(expr.Meta.Sample)
		fprintf(out, '\t%s', expr.Meta.Sample{k});
	end
	fprintf(out, '\n');
	
	for g = 1:length(organism.Genes.Name)
		fprintf(out, '%s', organism.Genes.Name{g});
		for s = 1:length(expr.Meta.Sample)
			fprintf(out, '\t%f', expr.Mean(g, s));
		end
		fprintf(out, '\n');
		progress.update(g / length(organism.Genes.Name));
	end
else
	error(['Exporting to tab delimited format is not supported for this ' ...
	       'type of data.']);
end

fclose(out);

