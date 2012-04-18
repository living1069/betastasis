function [] = xls_intersect(file_a, file_b)

[~, ~, raw] = xlsread(file_a);
a.ensembl = raw(2:end, 1);
a.hgnc = raw(2:end, 2);
a.entrez = raw(2:end, 3);
a.fold_change = str2double(raw(2:end, 4));

valid = false(length(a.ensembl), 1);
for k = 1:length(valid), valid(k) = ~all(isnan(a.ensembl{k})); end
a.ensembl = a.ensembl(valid);
a.hgnc = a.hgnc(valid);
a.entrez = a.entrez(valid);
a.fold_change = a.fold_change(valid);





[~, ~, raw] = xlsread(file_b);
b.ensembl = raw(2:end, 1);
b.entrez = raw(2:end, 3);
b.fold_change = str2double(raw(2:end, 4));

valid = false(length(b.ensembl), 1);
for k = 1:length(valid), valid(k) = ~all(isnan(b.ensembl{k})); end
b.ensembl = b.ensembl(valid);
b.entrez = b.entrez(valid);
b.fold_change = b.fold_change(valid);







if length(a.ensembl) ~= length(unique(a.ensembl))
	for k = 1:length(a.ensembl)
		if sum(strcmp(a.ensembl{k}, a.ensembl)) > 1
			fprintf('Not unique: %s\n', a.ensembl{k});
		end
	end
end
if length(b.ensembl) ~= length(unique(b.ensembl))
	for k = 1:length(b.ensembl)
		if sum(strcmp(b.ensembl{k}, b.ensembl)) > 1
			fprintf('Not unique: %s\n', b.ensembl{k});
		end
	end
end

[~, pos] = ismember(a.ensembl, b.ensembl);

fprintf('\n\n\nEnsembl\tHGNC\tEntrez\tlogFC (A)\tlogFC (B)\n');
for k = find(pos ~= 0)'
	fprintf('%s\t%s\t%s\t%.2f\t%.2f\n', a.ensembl{k}, a.hgnc{k}, ...
		num2str(a.entrez{k}), a.fold_change(k), b.fold_change(pos(k)));
end

