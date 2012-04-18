
% Author: Matti Annala <matti.annala@tut.fi>

function [] = mirna_similarity(trim_len)

global organism;
mirnas = organism.miRNA;

M = length(mirnas.Name);

valid = find(cellfun(@length, mirnas.Sequence) >= trim_len);
fprintf('Found %d / %d miRNA shorter than %d bases.\n', ...
	M - length(valid), M, trim_len);
	
mirna_trimmed = mirnas.Sequence(valid);
M = length(valid);
for m = 1:M, mirna_trimmed{m} = mirna_trimmed{m}(1:trim_len); end
	
seq_dist = zeros(M, M);
for m = 1:M
	for n = 1:m-1
		seq_dist(m, n) = sum(mirna_trimmed{n} ~= mirna_trimmed{m});
	end
end

seq_dist_vec = squareform(seq_dist);

for mm = 0:2
	[i, j] = find(tril(squareform(seq_dist_vec == mm)));
	
	fprintf('Found %d miRNA pairs with %d mismatches in %d first bases:\n', ...
		length(i), mm, trim_len);
	for k = 1:length(i)
		fprintf('- %s vs %s\n', mirnas.Name{valid(i(k))}, ...
			mirnas.Name{valid(j(k))});
	end

	fprintf('\n\n\n');
end

