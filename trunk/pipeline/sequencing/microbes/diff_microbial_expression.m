
function [] = diff_microbial_expression(microbes, test, ref, varargin)

if islogical(test), test = find(test); test = test(:)'; end
if islogical(ref), ref = find(ref); ref = ref(:)'; end

S = size(microbes.expr, 2);

%p = nan(length(microbes.name), 1);
%for v = 1:length(microbes.name)
%	p(v) = ranksum(microbes.expr(v, test), microbes.expr(v, ref));
%end
[~, p] = ttest2(microbes.expr(:, test)', microbes.expr(:, ref)');
p = p(:);

sample_order = [test, ref];

avg_viral_expr = nan(length(microbes.name), 2);
avg_viral_expr = [mean(microbes.expr(:, test), 2) mean(microbes.expr(:, ref), 2)];

[~, order] = sort(p, 'ascend');

fid = fopen('diff_microbial_expression.txt', 'W');

fprintf(fid, 'Microbe\tRanksum p\tAvg test expr\tAvg ref expr');
for s = sample_order
	if isfield(microbes.meta, 'sample_type')
		fprintf(fid, '\t%s (%s)', microbes.meta.sample_type{s}, ...
			microbes.meta.sample_id{s});
	else
		fprintf(fid, '\t%s', microbes.meta.sample_id{s});
	end
end
fprintf(fid, '\n');

for k = order'
	if p(k) > 0.05, continue, end
	
	fprintf(fid, '%s\t%.3f\t%.1f\t%.1f', ...
		regexprep(microbes.name{k}, '^gi\|.+\| ', ''), p(k), ...
		avg_viral_expr(k, 1), avg_viral_expr(k, 2));
	fprintf(fid, '\t%.1f', microbes.expr(k, sample_order));
	fprintf(fid, '\n');
end

fclose(fid);

