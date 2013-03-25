
function expr = pset_summary_median_polish(raw)

S = size(raw.mean, 2);

field = '';
if isfield(raw.rows, 'gene'), field = 'gene'; end
features = eval(sprintf('raw.rows.%s', field));

% Reorder the probes so that probes for each feature are placed contiguously.
[features, order] = sort(features);
raw = filter_rows(raw, order);

run_ends = [find(~strcmp(features(1:end-1), features(2:end)));length(features)];
run_lengths = diff([0; run_ends]);

expr = raw;
eval(sprintf('expr.rows.%s = features(run_ends);', field));
expr.rows.num_probes = run_lengths;

probe_indices = zeros(size(raw.mean, 1), 1);
pos = 1;
for r = 1:length(run_ends)
	probe_indices(pos:pos+run_lengths(r)-1) = 0:run_lengths(r)-1;
	pos = pos + run_lengths(r);
end

expr.mean = rmasummary(probe_indices, raw.mean, 'Output', 'natural');

