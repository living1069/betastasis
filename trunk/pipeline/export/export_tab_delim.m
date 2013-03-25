function [] = export_tab_delim(ds, file_prefix)

[file_prefix '.txt']
fid = fopen([file_prefix '.txt'], 'W')

if rx(ds.meta.type, '(gene|mirna|protein) expression')
	
	fprintf(fid, '\t%s', ds.meta.sample_id{:});
	fprintf(fid, '\n');
	
	if rx(ds.meta.type, 'gene expression')
		features = ds.rows.gene_symbol;
	elseif rx(ds.meta.type, 'pre-mirna expression')
		features = ds.rows.premirna_symbol;
	elseif rx(ds.meta.type, '^mirna expression')
		features = ds.rows.mirna_symbol;
	elseif rx(ds.meta.type, '^protein expression')
		features = ds.rows.protein;
	end
	
	for g = 1:length(features)
		fprintf(fid, '%s', features{g});
		fprintf(fid, '\t%.2f', ds.mean(g, :));
		fprintf(fid, '\n');
	end
else
	fprintf('Dataset type is not supported. Only exporting clinical data...\n');
end

fclose(fid);

if isfield(ds, 'meta')
	export_meta(ds.meta, [file_prefix '_clinical.txt']);
else
	fprintf('Clinical data not found.\n');
end









function [] = export_meta(meta, filename)

fid = fopen(filename, 'W');

fields = fieldnames(meta);
discard = strcmpi('type', fields);
fields = fields(~discard);

up_fields = upper(fields);
fprintf(fid, '%s\t', up_fields{1:end-1});
fprintf(fid, '%s\n', up_fields{end});

data = struct2cell(meta);
data = data(~discard);

for s = 1:length(meta.sample_id)
	for f = 1:length(fields)
		if iscellstr(data{f})
			fprintf(fid, '%s', data{f}{s});
		elseif isfloat(data{f})
			fprintf(fid, '%.2f', data{f}(s));
		elseif isinteger(data{f})
			fprintf(fid, '%d', data{f}(s));
		end
		if f < length(fields), fprintf(fid, '\t'); else fprintf(fid, '\n'); end
	end
end

fclose(fid);

