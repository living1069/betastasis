function tsp = find_tsp(expr, ref, features)

global organism;

genes = find(~isnan(sum([sum(expr.Mean, 2) sum(ref.Mean, 2)], 2)));
expr = expr.Mean(genes, :);
ref = ref.Mean(genes, :);

N = size(expr, 1);

pair = zeros(N * (N - 1) / 2, 2);
score = zeros(N * (N - 1) / 2, 1);

expr_rank = zeros(size(expr)); 
ref_rank = zeros(size(ref));

for k = 1:size(expr, 2)
	[~, order] = sort(expr(:, k), 'descend');
	expr_rank(order, k) = 1:length(order);
end

for k = 1:size(ref, 2)
	[~, order] = sort(ref(:, k), 'descend');
	ref_rank(order, k) = 1:length(order);
end

k = 1;
for a = 1:size(N)
	for b = a+1:size(N)
		p_expr = sum(expr_rank(a, :) < expr_rank(b, :)) / size(expr_rank, 2);
		ref_expr = sum(ref_rank(a, :) < ref_rank(b, :)) / size(ref_rank, 2);
		pair(k, :) = [a b];
		score(k) = abs(p_expr - ref_expr);
		k = k + 1;
	end
end

[~, order] = sort(score, 'descend');




% Find tied scores and reorder based on average rank difference between the
% two classes.
run_ends = [ find(score(1:end-1) ~= score(2:end)); length(score) ];
run_lengths = diff([0; run_ends]);

pos = 1;
for r = 1:length(run_lengths)
	if run_lengths(r) > 1
		sec_scores = zeros(run_lengths(r), 1);
		for k = 1:length(sec_scores)
			idx = pos + k - 1;
			expr_sec_score = ...
				mean(expr_rank(pair(idx, 1), :) - expr_rank(pair(idx, 2), :));
			ref_sec_score = ...
				mean(ref_rank(pair(idx, 1), :) - ref_rank(pair(idx, 2), :));
			sec_score(k) = abs(expr_sec_score - ref_sec_score);
		end
		
		[~, order] = sort(sec_score, 'descend');
		ref_rank = 
	end
	
	pos = pos + run_lengths(r);
end







