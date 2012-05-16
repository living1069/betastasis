function [] = export_graph_json(out_file, connectivity, labels, groups)

% Discard nodes with no connections.
keep = true(1, length(labels));
for k = 1:length(labels)
	if all(connectivity(k, :) == 0) && all(connectivity(:, k) == 0)
		keep(k) = false;
	end
end

labels = labels(keep);
connectivity = connectivity(keep, keep);
groups = groups(keep, keep);

fid = fopen(out_file, 'W')
fprintf(fid, '{"nodes":[');
for k = 1:length(labels)
	fprintf(fid, '{"name":"%s"}', labels{k});
	if k ~= length(labels), fprintf(fid, ','); end
end
fprintf(fid, '],\n"links":[');
[row, col] = find(connectivity > 0);
for k = 1:length(row)
	fprintf(fid, '{"source":%d, "target":%d, "value":%.1f, "type":%d}', ...
		row(k)-1, col(k)-1, connectivity(row(k), col(k)), ...
		groups(row(k), col(k)));
	if k ~= length(row), fprintf(fid, ','); end
end
fprintf(fid, ']}\n');
fclose(fid);

