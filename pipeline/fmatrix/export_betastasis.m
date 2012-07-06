function [] = export_betastasis(fmatrix, path)

global organism;

featidx = find(rx(fmatrix.features, 'N:EXPR:'));
tmp = fmatrix.data(featidx, :);
tmp(isnan(tmp)) = 0;   % FIXME: Maybe use something else?
fmatrix.data(featidx, :) = tmp;

expr_names = {};
tokens = regexpi(fmatrix.features, 'N:EXPR:(.+?)(:|$)', 'tokens');
for f = 1:length(fmatrix.features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	if rx(name, '/'), continue, end
	
	expr_names{end+1, 1} = name;
	
	dir = [path '/expr/' lower(name(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' name '.json'], 'data', fmatrix.data(f, :));
end

if length(expr_names) > 0
	export_json([path '/expr/features.json'], ...
		'features', expr_names);
end







featidx = find(rx(fmatrix.features, 'N:DIFF_EXPR:'));
tmp = fmatrix.data(featidx, :);
tmp(isnan(tmp)) = -100;   % FIXME: Maybe use something else?
fmatrix.data(featidx, :) = tmp;

expr_names = {};
tokens = regexpi(fmatrix.features, 'N:DIFF_EXPR:(.+?):', 'tokens');
for f = 1:length(fmatrix.features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	expr_names{end+1, 1} = name;
	
	dir = [path '/diff_expr/' lower(name(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' name '.json'], 'data', fmatrix.data(f, :));
end

if length(expr_names) > 0
	export_json([path '/diff_expr/features.json'], ...
		'features', expr_names);
end


	
	
	
	



featidx = find(rx(fmatrix.features, 'N:EXON_EXPR:'));
splice_genes = cell(length(featidx), 1);
splice_exons = cell(length(featidx), 1);

tokens = regexpi(fmatrix.features(featidx), 'N:EXON_EXPR:(.+)\[(.+)\]', 'tokens');
for f = 1:length(featidx)
	splice_genes{f} = tokens{f}{1}{1};
	splice_exons{f} = tokens{f}{1}{2};
end

tmp = fmatrix.data(featidx, :);
tmp(isnan(tmp)) = -100;   % FIXME: Maybe use something else?
fmatrix.data(featidx, :) = tmp;

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
		'data', fmatrix.data(idx, :));
end

if length(splice_exons) > 0
	export_json([path '/exon_expr/features.json'], ...
		'features', unique(splice_genes));
end
























featidx = find(rx(fmatrix.features, 'N:ESPL:'));
splice_genes = cell(length(featidx), 1);
splice_exons = cell(length(featidx), 1);

tokens = regexpi(fmatrix.features(featidx), 'N:ESPL:(.+)\[(.+)\]', 'tokens');
for f = 1:length(featidx)
	splice_genes{f} = tokens{f}{1}{1};
	splice_exons{f} = tokens{f}{1}{2};
end

tmp = fmatrix.data(featidx, :);
tmp(isnan(tmp)) = 0;   % FIXME: Maybe use something else?
fmatrix.data(featidx, :) = tmp;

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
		'data', fmatrix.data(idx, :));
end

if length(splice_exons) > 0
	export_json([path '/exon_splice/features.json'], ...
		'features', unique(splice_genes));
end

	
	

	
cna_names = {};
tokens = regexpi(fmatrix.features, 'N:CNA:(.+?):', 'tokens');
for f = 1:length(fmatrix.features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	cna_names{end+1, 1} = name;
	
	dir = [path '/cna/' lower(name(1))];
	[~, ~] = mkdir(dir);
	
	fmatrix.data(f, isnan(fmatrix.data(f, :))) = -100;
	export_json([dir '/' name '.json'], 'data', fmatrix.data(f, :));
end

if length(cna_names) > 0
	export_json([path '/cna/features.json'], 'features', cna_names);
end







mut_names = {};
tokens = regexpi(fmatrix.features, 'N:MUT:(.+?):', 'tokens');
for f = 1:length(fmatrix.features)
	if isempty(tokens{f}), continue, end
	
	name = tokens{f}{1}{1};
	mut_names{end+1, 1} = name;
	
	fmatrix.data(f, isnan(fmatrix.data(f, :))) = -100;
	
	dir = [path '/mutation/' lower(name(1))];
	[~, ~] = mkdir(dir);
	export_json([dir '/' name '.json'], 'data', fmatrix.data(f, :));
end

if length(mut_names) > 0
	export_json([path '/mutation/features.json'], 'features', mut_names);
end



	

	
	
% WRITE THE NEW-STYLE CLINICAL DATA
featidx = find(rx(fmatrix.features, '.:CLIN:'));
tmp = fmatrix.data(featidx, :);
tmp(isnan(tmp)) = -1;
fmatrix.data(featidx, :) = tmp;

clin_args = {'sample_id', fmatrix.samples};

tokens = regexpi(fmatrix.features, 'N:CLIN:(.+)', 'tokens');
for f = 1:length(fmatrix.features)
	if isempty(tokens{f}), continue, end
	
	clin_args{end+1} = tokens{f}{1}{1};
	clin_args{end+1} = fmatrix.data(f, :);
end

tokens = regexpi(fmatrix.features, 'C:CLIN:([^:]+):(.+)', 'tokens');
for f = 1:length(fmatrix.features)
	if isempty(tokens{f}), continue, end
		
	cat_names = textscan(tokens{f}{1}{2}, '%s', 'Delimiter', ';');
	cat_names = cat_names{1};
	
	clin_args{end+1} = tokens{f}{1}{1};
	clin_args{end+1} = cat_names(fmatrix.data(f, :));
end

export_json([path '/clinical.json'], clin_args{:});

	
