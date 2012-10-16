function hugo = ensembl_to_hugo(gene_ids, gtf_file)

data = readtable(gtf_file, 'Ignore', [1:8 10:100], 'Header', false);
tokens = rx_capture(data{1}, 'gene_id "(.*?)";.*; gene_name "(.*?)"');

[~, idx] = unique(tokens(:, 1));
tokens = tokens(idx, :);

hugo = repmat({''}, length(gene_ids), 1);
[~, pos] = ismember(gene_ids, tokens(:, 1));

hugo(pos ~= 0) = tokens(pos(pos ~= 0), 2);

