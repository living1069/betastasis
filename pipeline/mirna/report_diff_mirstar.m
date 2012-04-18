function [] = report_diff_mirstar(test, ref, varargin)

global organism;
mirnas = organism.miRNA;

min_expr = -Inf;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MinExpr')
		min_expr = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isfield(test, 'miRNA'), test_expr = test.miRNA;
else test_expr = test.Mean; end
	
if isfield(ref, 'miRNA'), ref_expr = ref.miRNA;
else ref_expr = ref.Mean; end

M = length(mirnas.Name);

% First we need to find all functional and degraded miRNA pairs.
mirstar_pair = zeros(0, 2);
for m = 1:M
	star = mirna_idx([mirnas.Name{m} '*']);
	if isnan(star), continue, end
	mirstar_pair(end+1, :) = [m star];
end

valid = ~all(isnan(test_expr(mirstar_pair(:, 1), :)), 2) & ...
	~all(isnan(test_expr(mirstar_pair(:, 2), :)), 2) & ...
	~all(isnan(ref_expr(mirstar_pair(:, 1), :)), 2) & ...
	~all(isnan(ref_expr(mirstar_pair(:, 2), :)), 2);
	
valid = valid & ...
	(median(test_expr(mirstar_pair(:, 1), :), 2) >= min_expr | ...
	median(test_expr(mirstar_pair(:, 2), :), 2) >= min_expr) & ...
	(median(ref_expr(mirstar_pair(:, 1), :), 2) >= min_expr | ...
	median(ref_expr(mirstar_pair(:, 2), :), 2) >= min_expr);

mirstar_pair = mirstar_pair(valid, :);

% We add a small tweak factor to all the expression values to prevent
% saturation effects like infinite ratios.
test_expr = test_expr + 0.5;
ref_expr = ref_expr + 0.5;

% Then we calculate the miR/miR* ratio in each sample, and look for the
% miRNA with the highest differential ratio between test and reference samples.
test_ratios = nanmedian(test_expr(mirstar_pair(:, 1), :) ./ ...
	test_expr(mirstar_pair(:, 2), :), 2);
ref_ratios = nanmedian(ref_expr(mirstar_pair(:, 1), :) ./ ...
	ref_expr(mirstar_pair(:, 2), :), 2);
	
diff_ratios = log2(test_ratios ./ ref_ratios);

fid = fopen('diff_mirstar.txt', 'W');
[~, order] = sort(abs(diff_ratios), 'descend');
for m = order'
	fprintf(fid, '%s\t%f\n', mirnas.Name{mirstar_pair(m, 1)}, diff_ratios(m));
end
fclose(fid);

