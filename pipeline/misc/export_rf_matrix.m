function [] = export_rf_matrix(fmatrix, path)

% Impute numeric values from nearest neighbor columns. For categorical
% features, we change unknown values to -1.
numeric = true(length(fmatrix.Features), 1);
for k = 1:length(numeric)
	if fmatrix.Features{k}(1) ~= 'N', numeric(k) = false; end
end

fmatrix.Data(numeric, :) = knnimpute(fmatrix.Data(numeric, :), 3);
categoric = fmatrix.Data(~numeric, :);
categoric(isnan(categoric)) = -1;
fmatrix.Data(~numeric, :) = categoric;

out_tmp = ptemp;
dlmwrite(out_tmp, fmatrix.Data', '\t');

tmp_fid = fopen(out_tmp);
final_fid = fopen(path, 'W');

for f = 1:length(fmatrix.Features)
	fprintf(final_fid, '\t%s', fmatrix.Features{f});
end
fprintf(final_fid, '\n');

for s = 1:length(fmatrix.Samples)
	line = fgetl(tmp_fid);
	if line == -1, error 'Data file ended abruptly.', end
	
	fprintf(final_fid, '%s\t%s\n', fmatrix.Samples{s}, line);
end

fclose(final_fid);
fclose(tmp_fid);

safe_delete(out_tmp);

