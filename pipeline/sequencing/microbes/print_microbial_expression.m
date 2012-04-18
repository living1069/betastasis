
function [] = print_microbial_expression(microbes, varargin)

min_reads = 10;
min_unique_read_ratio = 0.2;
min_unique_reads = 50;

S = size(microbes.expr, 2);

fid = fopen('microbial_expression.txt', 'W');

if isfield(microbes.meta, 'sample_type')
	[~, order] = sort_nat(microbes.meta.sample_type);
else
	order = (1:S)';
end

for s = order'
	[~, order] = sort(microbes.expr(:, s), 'descend');
	
	fprintf(fid, '\nSAMPLE %s:\n', microbes.meta.sample_id{s});
	for k = order'
		num_reads = microbes.expr(k, s);
		num_unique_reads = length(unique(microbes.read_seqs{k, s}));
		
		if num_reads <= min_reads, continue, end
		
		if num_unique_reads <= min_unique_reads && ...
			num_unique_reads / num_reads <= min_unique_read_ratio
			continue;
		end
		
		fprintf(fid, '%s: %d reads (%d unique)\n', microbes.name{k}, ...
			num_reads, num_unique_reads);
	end
end

fclose(fid);

