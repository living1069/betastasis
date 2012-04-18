function [] = report_isomirs(isomirs, report, varargin)

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

groups = [];

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Groups')
		groups = varargin{k+1};
		if length(groups) > 2, error 'Too many groups.'; end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

M = length(mirnas.Name);
S = size(isomirs.Isoforms, 2);

fid = fopen(report, 'W');

for s = 1:S
	fprintf(fid, '\t%s', isomirs.Meta.Sample.Filename{s});
end
fprintf(fid, '\tTotal\n');

if ~isempty(groups)
	divergence = nan(1, M);
	
	for m = 1:M
		isomir_expr = zeros(length(isomirs.Isoforms{m, 1}.Count), S);
		for s = 1:S
			isomir_expr(:, s) = isomirs.Isoforms{m, s}.Count;
		end
		
		% For calculating the divergences, we require that the miRNA must have
		% an adequate number of reads in both groups.
		valid = find(min(mean(isomir_expr(:, groups{1}), 2), ...
			mean(isomir_expr(:, groups{2}), 2)) > 10);
		isomir_expr = isomir_expr(valid, :);
		
		if numel(isomir_expr) == 0, continue, end
		
		g1_median = median(isomir_expr(:, groups{1}), 2);
		g2_median = median(isomir_expr(:, groups{2}), 2);
		
		% L1 divergence between the probability distributions.
		% FIXME: Should we take the divergence over the cdf instead.
		divergence(m) = sum(abs(g1_median / sum(g1_median) - ...
			g2_median / sum(g2_median)));
	end
	
	divergence(isnan(divergence)) = 0;
	[~, mirna_order] = sort(divergence, 'descend');
	
	divergence(mirna_order(1:5))
	divergence(mirna_order(end-4:end))

else
	mirna_order = 1:M;
end

for m = mirna_order
	if length(isomirs.Isoforms{m, 1}.Sequence) ~= ...
		length(isomirs.Isoforms{m, 1}.Count)
		continue;
	end
	
	isomir_expr = zeros(length(isomirs.Isoforms{m, 1}.Count), S);
	for s = 1:S
		isomir_expr(:, s) = isomirs.Isoforms{m, s}.Count;
	end
	
	sum_expr = sum(isomir_expr, 2);
	
	valid = find(sum_expr / S > 10);
	[~, order] = sort(sum_expr(valid), 'descend');
	order = valid(order);
	
	if isempty(order), continue, end
	
	fprintf(fid, '%s\n', mirnas.Name{m});
	
	for k = order'
		fprintf(fid, '%s', isomirs.Isoforms{m, s}.Sequence{k});
		fprintf(fid, '\t%d', isomir_expr(k, :));
		fprintf(fid, '\t%d', sum_expr(k));
		
		if strcmp(isomirs.Isoforms{m, s}.Sequence{k}, mirnas.Sequence{m})
			fprintf(fid, '\t***');
		end
		
		fprintf(fid, '\n');
	end
	
	% Calculate uridylation levels.
	uridylated = zeros(1, S);
	for k = 1:size(isomir_expr, 1)
		seq = isomirs.Isoforms{m, 1}.Sequence{k};
		if seq(end) == 'T'
			uridylated = uridylated + isomir_expr(k, :);
		end
	end
	
	uridylated = uridylated ./ sum(isomir_expr, 1);
	fprintf(fid, 'Uridylation level');
	fprintf(fid, '\t%.1f%%', uridylated * 100);
	fprintf(fid, '\n');
	
	% Calculate adenylation levels.
	adenylated = zeros(1, S);
	for k = 1:size(isomir_expr, 1)
		seq = isomirs.Isoforms{m, 1}.Sequence{k};
		if seq(end) == 'A'
			adenylated = adenylated + isomir_expr(k, :);
		end
	end
	
	adenylated = adenylated ./ sum(isomir_expr, 1);
	fprintf(fid, 'Adenylation level');
	fprintf(fid, '\t%.1f%%', adenylated * 100);
	fprintf(fid, '\n');

	
	fprintf(fid, '\n');
end

fclose(fid);

