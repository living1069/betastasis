function [] = export_betastasis(fmatrix, path)

global organism;

featidx = find(rx(fmatrix.Features, 'N:EXPR:'));
tmp = fmatrix.Data(featidx, :);
tmp(isnan(tmp)) = 0;   % FIXME: Maybe use something else?
fmatrix.Data(featidx, :) = tmp;

expr_names = {};
tokens = regexpi(fmatrix.Features, 'N:EXPR:(.+?):', 'tokens');
for f = 1:length(fmatrix.Features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	expr_names{end+1, 1} = name;
	
	dir = [path '/expr/' lower(name(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' name '.json'], 'data', fmatrix.Data(f, :));
end

if length(expr_names) > 0
	export_json([path '/expr/features.json'], ...
		'features', expr_names);
end







featidx = find(rx(fmatrix.Features, 'N:DIFF_EXPR:'));
tmp = fmatrix.Data(featidx, :);
tmp(isnan(tmp)) = -100;   % FIXME: Maybe use something else?
fmatrix.Data(featidx, :) = tmp;

expr_names = {};
tokens = regexpi(fmatrix.Features, 'N:DIFF_EXPR:(.+?):', 'tokens');
for f = 1:length(fmatrix.Features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	expr_names{end+1, 1} = name;
	
	dir = [path '/diff_expr/' lower(name(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' name '.json'], 'data', fmatrix.Data(f, :));
end

if length(expr_names) > 0
	export_json([path '/diff_expr/features.json'], ...
		'features', expr_names);
end


	
	
	
	



featidx = find(rx(fmatrix.Features, 'N:EXON_EXPR:'));
splice_genes = cell(length(featidx), 1);
splice_exons = cell(length(featidx), 1);

tokens = regexpi(fmatrix.Features(featidx), 'N:EXON_EXPR:(.+)\[(.+)\]', 'tokens');
for f = 1:length(featidx)
	splice_genes{f} = tokens{f}{1}{1};
	splice_exons{f} = tokens{f}{1}{2};
end

tmp = fmatrix.Data(featidx, :);
tmp(isnan(tmp)) = 0;   % FIXME: Maybe use something else?
fmatrix.Data(featidx, :) = tmp;

done = false(length(featidx), 1);
for f = 1:length(featidx)
	if done(f), continue, end
	
	gene = splice_genes{f};
	ex = find(strcmp(gene, splice_genes));
	
	done(ex) = true;
	
	% Perform natural sorting on the exons (i.e. 2 comes before 10)
	[~, order] = sort_nat(splice_exons(ex));
	ex = ex(order);
	
	idx = featidx(ex);
	
	dir = [path '/exon_expr/' lower(gene(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' gene '.json'], ...
		'exons', splice_exons(ex), ...
		'data', fmatrix.Data(idx, :));
end

if length(splice_exons) > 0
	export_json([path '/exon_expr/features.json'], ...
		'features', unique(splice_genes));
end
























featidx = find(rx(fmatrix.Features, 'N:ESPL:'));
splice_genes = cell(length(featidx), 1);
splice_exons = cell(length(featidx), 1);

tokens = regexpi(fmatrix.Features(featidx), 'N:ESPL:(.+)\[(.+)\]', 'tokens');
for f = 1:length(featidx)
	splice_genes{f} = tokens{f}{1}{1};
	splice_exons{f} = tokens{f}{1}{2};
end

tmp = fmatrix.Data(featidx, :);
tmp(isnan(tmp)) = 0;   % FIXME: Maybe use something else?
fmatrix.Data(featidx, :) = tmp;

done = false(length(featidx), 1);
for f = 1:length(featidx)
	if done(f), continue, end
	
	gene = splice_genes{f};
	ex = find(strcmp(gene, splice_genes));
	
	done(ex) = true;
	
	% Perform natural sorting on the exons (i.e. 2 comes before 10)
	[~, order] = sort_nat(splice_exons(ex));
	ex = ex(order);
	
	idx = featidx(ex);
	
	dir = [path '/exon_splice/' lower(gene(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' gene '.json'], ...
		'exons', splice_exons(ex), ...
		'data', fmatrix.Data(idx, :));
end

if length(splice_exons) > 0
	export_json([path '/exon_splice/features.json'], ...
		'features', unique(splice_genes));
end

	
	

	
cna_names = {};
tokens = regexpi(fmatrix.Features, 'N:CNA:(.+?):', 'tokens');
for f = 1:length(fmatrix.Features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	cna_names{end+1, 1} = name;
	
	dir = [path '/cna/' lower(name(1))];
	[~, ~] = mkdir(dir);
	
	fmatrix.Data(f, isnan(fmatrix.Data(f, :))) = -100;
	export_json([dir '/' name '.json'], 'data', fmatrix.Data(f, :));
end

if length(cna_names) > 0
	export_json([path '/cna/features.json'], 'features', cna_names);
end







mut_names = {};
tokens = regexpi(fmatrix.Features, 'N:MUT:(.+?):', 'tokens');
for f = 1:length(fmatrix.Features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	mut_names{end+1, 1} = name;
	
	fmatrix.Data(f, isnan(fmatrix.Data(f, :))) = -100;
	
	dir = [path '/mutation/' lower(name(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' name '.json'], 'data', fmatrix.Data(f, :));
end

if length(mut_names) > 0
	export_json([path '/mutation/features.json'], 'features', mut_names);
end



	


	
dir = [path '/clinical'];
[~, ~] = mkdir(dir);

featidx = find(rx(fmatrix.Features, '.:CLIN:'));
tmp = fmatrix.Data(featidx, :);
tmp(isnan(tmp)) = -1;
fmatrix.Data(featidx, :) = tmp;

clinical_features = {};
tokens = regexpi(fmatrix.Features, 'N:CLIN:(.+)', 'tokens');
for f = 1:length(fmatrix.Features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	clinical_features{end+1, 1} = name;
	
	export_json([dir '/' name '.json'], 'data', fmatrix.Data(f, :));
end

tokens = regexpi(fmatrix.Features, 'C:CLIN:([^:]+):(.+)', 'tokens');
for f = 1:length(fmatrix.Features)
	if isempty(tokens{f}), continue, end
		
	cat_names = textscan(tokens{f}{1}{2}, '%s', 'Delimiter', ',');
	cat_names = cat_names{1}
	
	name = tokens{f}{1}{1};
	clinical_features{end+1, 1} = name;
	
	export_json([dir '/' name '.json'], 'data', cat_names(fmatrix.Data(f, :)));
end


if ~isempty(clinical_features)
	export_json([path '/clinical/features.json'], ...
		'features', clinical_features);
end

export_json([path '/clinical/sample_id.json'], ...
	'data', fmatrix.Samples);

