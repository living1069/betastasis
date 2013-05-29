
function expr = pset_summary_median(raw)

S = size(raw.mean, 2);

field = '';
if isfield(raw.rows, 'gene'), field = 'gene'; end
if isfield(raw.rows, 'target'), field = 'target'; end
features = eval(sprintf('raw.rows.%s', field));

[features, order] = sort(features);
raw = filter_rows(raw, order);

run_ends = [find(~strcmp(features(1:end-1), features(2:end)));length(features)];
run_lengths = diff([0; run_ends]);

expr = raw;
eval(sprintf('expr.rows.%s = features(run_ends);', field));
expr.rows.num_probes = run_lengths;
expr.mean = nan(length(run_ends), size(raw.mean, 2));

pos = 1;
for r = 1:length(run_ends)
	expr.mean(r, :) = nanmedian(raw.mean(pos:pos+run_lengths(r)-1, :), 1);
	pos = pos + run_lengths(r);
end


