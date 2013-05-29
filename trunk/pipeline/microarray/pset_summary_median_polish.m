
function expr = pset_summary_median_polish(raw)

S = size(raw.mean, 2);

field = '';
if isfield(raw.rows, 'name'), field = 'name'; end
if isfield(raw.rows, 'gene'), field = 'gene'; end
targets = eval(sprintf('raw.rows.%s', field));

% Throw away all rows that do not correspond to any targets.
valid = ~strcmp(targets, '');
raw = filter_rows(raw, valid);
targets = targets(valid);

% Reorder the probes so that probes for each feature are placed contiguously.
[targets, order] = sort(targets);
raw = filter_rows(raw, order);

run_ends = [find(~strcmp(targets(1:end-1), targets(2:end)));length(targets)];
run_lengths = diff([0; run_ends]);

expr = raw;
eval(sprintf('expr.rows.%s = targets(run_ends);', field));
expr.rows.num_probes = run_lengths;

probe_indices = zeros(size(raw.mean, 1), 1);
pos = 1;
for r = 1:length(run_ends)
	probe_indices(pos:pos+run_lengths(r)-1) = 0:run_lengths(r)-1;
	pos = pos + run_lengths(r);
end

expr.mean = rmasummary(probe_indices, raw.mean, 'Output', 'natural');

