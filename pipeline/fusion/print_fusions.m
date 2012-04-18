function fusions_ordered = print_fusions(rearrangements, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;
chromosomes = organism.Chromosomes;

ranking = 'chisq';
gene_expr = [];
exon_expr = [];
groups = {};

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Ranking')
		ranking = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if regexpi(varargin{k}, 'gene.*expr')
		gene_expr = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if regexpi(varargin{k}, 'exon.*expr')
		exon_expr = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if strcmpi(varargin{k}, 'Groups')
		groups = varargin{k+1};
		for g = 1:2:length(groups)
			if ~ischar(groups{g}), error 'Invalid group format.'; end
			if islogical(groups{g+1})
				groups{g+1} = find(groups{g+1});
			elseif ~isnumeric(groups{g+1})
				error 'Invalid group format.';
			end
		end
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);


fusions = rearrangements.fusions;
S = length(fusions);

joint_fusions = filter_fusions(rearrangements, varargin{:});



% Check if the given gene expression data matches with our fusion data.
if ~isempty(gene_expr)
	[found, idx] = ismember(rearrangements.meta.sample_id, ...
		gene_expr.meta.sample_id);
	if any(~found) || length(idx) ~= S
		error 'Expression data does not match with rearrangement samples.';
	end
	gene_expr = filter_query(gene_expr, idx);
end

if ~isempty(exon_expr)
	[found, idx] = ismember(rearrangements.meta.sample_id, ...
		exon_expr.meta.sample_id);
	if any(~found) || length(idx) ~= S
		error 'Expression data does not match with rearrangement samples.';
	end
	exon_expr = filter_query(exon_expr, idx);
end




if ~isempty(gene_expr) && ~isempty(exon_expr)
	error 'Combining gene and exon expression data is not supported yet.';
end

% Calculate differential expression if we have exon expression data available.
if ~isempty(exon_expr)
	
	% 5' of 5' gene, 3' of 5' gene, 5' of 3' gene, 3' of 3' gene
	diff_expr = nan(length(joint_fusions), 4);
	
	for f = 1:length(joint_fusions)
		gfusion = joint_fusions{f};
		
		% Sort the gene exons based on their order within the gene.
		exons_5p_gene = find(exons.Gene == gfusion.Genes(1));
		exons_3p_gene = find(exons.Gene == gfusion.Genes(2));
		
		[~, order] = sort_nat(exons.ID(exons_5p_gene));
		exons_5p_gene = exons_5p_gene(order)
		
		[~, order] = sort_nat(exons.ID(exons_3p_gene));
		exons_3p_gene = exons_3p_gene(order)

		neg_samples = true(1, S);
		for p = 1:size(gfusion.Exons, 1)
			neg_samples(gfusion.ReadSamples{p}) = false;
		end
		neg_samples
		
		% Convert the fusion exon IDs to indices against the vectors
		% exons_*p_gene, that have been sorted by chromosomal exon position.
		fus_exons = nan(size(gfusion.Exons));
		for p = 1:size(gfusion.Exons, 1)
			fus_exons(p, 1) = find(exons_5p_gene == gfusion.Exons(p, 1));
			fus_exons(p, 2) = find(exons_3p_gene == gfusion.Exons(p, 2));
		end
		fus_exons
		
		% Figure out the per-sample fusion boundaries. If any one sample has
		% more than one breakpoints, then we select maximal 5' and 3' regions
		% that do not overlap with any found breakpoints.
		% [5p last 5p exon, 5p first 3p exon, 3p last 5p exon, 3p first 3p exon]
		breaks = nan(S, 4);
		for p = 1:size(fus_exons, 1)
			s = unique(gfusion.ReadSamples{p});
			breaks(s, 1) = min(breaks(s, 1), fus_exons(p, 1));
			breaks(s, 2) = max(breaks(s, 2), fus_exons(p, 1));
			breaks(s, 3) = min(breaks(s, 3), fus_exons(p, 2));
			breaks(s, 4) = max(breaks(s, 4), fus_exons(p, 2));
		end
		breaks
		
		neg_exon_expr_5p = mean(exon_expr.Mean(exons_5p_gene, neg_samples),2);
		neg_exon_expr_3p = mean(exon_expr.Mean(exons_3p_gene, neg_samples),2);
		
		% Calculate per-sample differential expression for the 4 exon groups.
		per_sample_diff_expr = zeros(S, 4);
		for s = find(~neg_samples)
			per_sample_diff_expr(s, 1) = ...
				sum(exon_expr.Mean(exons_5p_gene(1:breaks(s, 1)), s)) / ...
				sum(neg_exon_expr_5p(1:breaks(s, 1)));
			per_sample_diff_expr(s, 2) = ...
				sum(exon_expr.Mean(exons_5p_gene(breaks(s, 2):end), s)) / ...
				sum(neg_exon_expr_5p(breaks(s, 2):end));
			per_sample_diff_expr(s, 3) = ...
				sum(exon_expr.Mean(exons_3p_gene(1:breaks(s, 3)), s)) / ...
				sum(neg_exon_expr_3p(1:breaks(s, 3)));
			per_sample_diff_expr(s, 4) = ...
				sum(exon_expr.Mean(exons_3p_gene(breaks(s, 4):end), s)) / ...
				sum(neg_exon_expr_3p(breaks(s, 4):end));
		end
		
		per_sample_diff_expr
		
		% Finally, we calculate the overall differential expression for each
		% of the four exon groups.
		diff_expr(f, :) = log2(median(per_sample_diff_expr(~neg_samples), 1));
	end
end








if regexpi(ranking, 'total')
	
	fusion_reads = zeros(length(joint_fusions), 1);
	for f = 1:length(joint_fusions)
		fusion_reads(f) = sum(cellfun(@length, joint_fusions{f}.ReadSamples));
	end
	[~, order] = sort(fusion_reads, 'descend');
	

elseif regexpi(ranking, 'chi-?sq')
	
	chisq = nan(size(joint_fusions));
	for f = 1:length(joint_fusions)
		gfusion = joint_fusions{f};
		num_exon_pairs = size(gfusion.Exons, 1);
		
		sr = zeros(1, S);
		for p = 1:num_exon_pairs
			for s = 1:S, sr(s) = sr(s) + sum(gfusion.ReadSamples{p} == s); end
		end
		expected_reads = sum(sr) / S;
		chisq(f) = sum((sr - ones(1, S) * expected_reads).^2 / expected_reads);
	end
	
	[~, order] = sort(chisq, 'descend');
	
	
elseif regexpi(ranking, 'mmes')
	
	mmes_scores = zeros(length(joint_fusions), 1);
	
	for f = 1:length(joint_fusions)
		gfusion = joint_fusions{f};

		% Check the minimum anchor length.
		if ~isfield(gfusion, 'ReadSequences')
			error 'MMES ranking requires read sequence information.';
		end
		
		for p = 1:size(gfusion.Exons, 1)
			[seqs, idx] = unique(gfusion.ReadSequences{p});
			junction_offsets = gfusion.JunctionOffsets{p}(idx);
			
			for s = 1:length(seqs)
				left_len = junction_offsets(s);
				right_len = length(seqs{s}) - left_len;
				mmes_scores(f) = mmes_scores(f) + min(left_len, right_len);
			end
		end
	end
	
	[~, order] = sort(mmes_scores, 'descend');
	
	
elseif regexpi(ranking, 'gene.expr')
	
	diff_expr = nan(length(joint_fusions), 2);
	
	for f = 1:length(joint_fusions)
		gfusion = joint_fusions{f};
		
		samples_with_fusion = false(S, 1);
		for p = 1:size(gfusion.Exons, 1)
			samples_with_fusion(gfusion.ReadSamples{p}) = true;
		end
		
		for g = 1:2
			diff_expr(f, g) = log2(median(expr.Mean( ...
				gfusion.Genes(g), samples_with_fusion)+1)) - ...
				log2(median(expr.Mean( ...
				gfusion.Genes(g), ~samples_with_fusion)+1));
		end
	end
	
	[~, order] = sort(max(abs(diff_expr), [], 2), 'descend');

elseif regexpi(ranking, 'exon.expr')
	
	score = nanmax([abs(diff_expr(:, 1) - diff_expr(:, 2)), ...
		abs(diff_expr(:, 3) - diff_expr(:, 4))], [], 2);
	[~, order] = sort(score, 'descend');
	order = order(sum(isnan(score))+1:end);

else
	error 'Requested ranking method is not supported.';
end







% Print the full fusion sequences.
fid = fopen('fusion_reads.txt', 'W');

for k = 1:length(order)
	idx = order(k);
	gfusion = joint_fusions{idx};
	
	total_reads = sum(cellfun(@length, gfusion.ReadSamples));
	
	fprintf(fid, '- fusion of %s and %s (%d total reads):\n', ...
		genes.Name{gfusion.Genes(1)}, genes.Name{gfusion.Genes(2)}, ...
		total_reads);
	
	num_exon_pairs = size(gfusion.Exons, 1);
	for p = 1:num_exon_pairs
		left_exon = gfusion.Exons(p, 1);
		right_exon = gfusion.Exons(p, 2);
		left_exon_seq = upper(exons.Sequence{left_exon});
		right_exon_seq = upper(exons.Sequence{right_exon});
		
		fprintf(fid, '  * junction between %s[%s] and %s[%s]:\n', ...
			genes.Name{gfusion.Genes(1)}, exons.ID{gfusion.Exons(p, 1)}, ...
			genes.Name{gfusion.Genes(2)}, exons.ID{gfusion.Exons(p, 2)});
		
		for s = 1:length(gfusion.ReadSequences{p})
			left_len = gfusion.JunctionOffsets{p}(s);
			seq = gfusion.ReadSequences{p}{s};
				
			ref_seq = [left_exon_seq(end-left_len+1:end) ...
				right_exon_seq(1:length(seq)-left_len)];
				
			mismatches = (seq ~= ref_seq);
			seq(mismatches) = lower(seq(mismatches));
				
			fprintf(fid, '    * %d: %s|', gfusion.ReadSamples{p}(s), ...
				seq(1:left_len));
			fprintf(fid, '%s\n', seq(left_len+1:end));
		end
	end
end

fclose(fid);





% Construct a spreadsheet of the discovered fusions.
fid = fopen('fusion_spreadsheet.txt', 'W');
fprintf(fid, 'Fusion\t5'' locus\t3'' locus\t5'' exon(s)\t3'' exon(s)\tJunction type\tFrameshift\tTotal positive\t');

if ~isempty(groups)
	for k = 1:2:length(groups)
		if ~ischar(groups{k}), error 'Invalid group format.'; end
		fprintf(fid, '%s group avg reads\t', groups{k}); 
	end
end

if ~isempty(gene_expr), fprintf(fid, 'DE1\tDE2\t'); end
if ~isempty(exon_expr), fprintf(fid, '5p5p\t5p3p\t3p5p\t3p3p\t'); end
for s = 1:S
	if isfield(rearrangements.meta, 'sample_type')
		fprintf(fid, '%s (%s)\t', rearrangements.meta.sample_type{s}, ...
			rearrangements.meta.sample_id{s});
	else
		fprintf(fid, '%s\t', rearrangements.meta.sample_id{s});
	end
end
fprintf(fid, '\n');

fusions_ordered = cell(length(order), 1);

for k = 1:length(order)
	idx = order(k);
	gfusion = joint_fusions{idx};
	
	left_gene = gfusion.Genes(1);
	right_gene = gfusion.Genes(2);
	
	left_chr = genes.Chromosome(left_gene);
	right_chr = genes.Chromosome(right_gene);
	
	left_pos = genes.Position(left_gene, 1);
	right_pos = genes.Position(right_gene, 1);
	
	left_locus = '???';
	right_locus = '???';
	if ~isnan(left_chr)
		left_locus = sprintf('chr%s: %d kb', chromosomes.Name{left_chr}, ...
			round(left_pos / 1000));
	end
	
	if ~isnan(right_chr)
		right_locus = sprintf('chr%s: %d kb', chromosomes.Name{right_chr}, ...
			round(right_pos / 1000));
	end
		
	left_strand = genes.Strand(left_gene);
	right_strand = genes.Strand(right_gene);
	
	left_locus = [left_locus ' (' left_strand ')'];
	right_locus = [right_locus ' (' right_strand ')'];
	
	
	
	fusions_ordered{k} = [genes.Name{left_gene} ':' genes.Name{right_gene}];
	
	
	
	% Now we attempt to calculate the biological impact of the fusion.
	left_txs = genes.Transcripts(left_gene, ...
		1:genes.TranscriptCount(left_gene));
	right_txs = genes.Transcripts(right_gene, ...
		1:genes.TranscriptCount(right_gene));
	
	junction_types = {};
	frameshifts = [];
	
	for e = 1:size(gfusion.Exons, 1)
		
		left_exon_roles = {};
		right_exon_roles = {};
		left_exon_frames = [];
		right_exon_frames = [];
		
		for tx = left_txs
			cds = transcripts.CDS(tx, :);
			if all(isnan(cds))
				left_exon_roles{end+1, 1} = 'NC'; continue;
			end
			
			ex = find(transcripts.Exons{tx} == gfusion.Exons(e, 1));
			if isempty(ex), continue, end
				
			exon_pos = transcripts.ExonPos{tx};
			exon_pos = exon_pos(ex, :);
			
			if exon_pos(2) < cds(1)
				left_exon_roles{end+1, 1} = '5'' UTR';
			elseif exon_pos(2) >= cds(1) && exon_pos(2) <= cds(2)
				left_exon_roles{end+1, 1} = 'CDS';
				left_exon_frames(end+1) = mod(exon_pos(2) - cds(1) + 1, 3);
			elseif exon_pos(2) > cds(2)
				left_exon_roles{end+1, 1} = '3'' UTR';
			end
		end
		
		for tx = right_txs
			cds = transcripts.CDS(tx, :);
			if all(isnan(cds))
				right_exon_roles{end+1, 1} = 'NC'; continue;
			end
			
			ex = find(transcripts.Exons{tx} == gfusion.Exons(e, 2));
			if isempty(ex), continue, end
				
			exon_pos = transcripts.ExonPos{tx};
			exon_pos = exon_pos(ex, :);
			
			if exon_pos(2) < cds(1)
				right_exon_roles{end+1, 1} = '5'' UTR';
			elseif exon_pos(2) >= cds(1) && exon_pos(2) <= cds(2)
				right_exon_roles{end+1, 1} = 'CDS';
				right_exon_frames(end+1) = mod(exon_pos(1) - cds(1), 3);
			elseif exon_pos(2) > cds(2)
				right_exon_roles{end+1, 1} = '3'' UTR';
			end
		end
		
		
		% If only some transcripts of a gene are non-coding, we ignore them.
		left_exon_roles = unique(left_exon_roles);
		right_exon_roles = unique(right_exon_roles);
		
		if any(strcmpi('NC', left_exon_roles)) && ...
			length(left_exon_roles) > 1
			left_exon_roles = ...
				left_exon_roles(~strcmpi('NC', left_exon_roles));
		end
		
		if any(strcmpi('NC', right_exon_roles)) && ...
			length(right_exon_roles) > 1
			right_exon_roles = ...
				right_exon_roles(~strcmpi('NC', right_exon_roles));
		end
		
		left_exon_frames = unique(left_exon_frames);
		right_exon_frames = unique(right_exon_frames);
		
		if length(left_exon_frames) == 1 && length(right_exon_frames) == 1
			frameshifts(end+1, 1) = ...
				mod(left_exon_frames - right_exon_frames, 3);
		end
		
		for l = 1:length(left_exon_roles)
			for r = 1:length(right_exon_roles)
				junction_types{end+1, 1} = ...
					[left_exon_roles{l} ' to ' right_exon_roles{r}];
			end
		end
	end
	
	junction_types = unique(junction_types);
	frameshifts = unique(frameshifts);
	
	fprintf(fid, '%s:%s\t%s\t%s\t', ...
		genes.Name{left_gene}, genes.Name{right_gene}, left_locus, right_locus);
	
	left_exons = {};
	right_exons = {};
	
	num_exon_pairs = size(gfusion.Exons, 1);
	for p = 1:num_exon_pairs
		left_exons{end+1} = exons.ID{gfusion.Exons(p, 1)};
		right_exons{end+1} = exons.ID{gfusion.Exons(p, 2)};
	end
	
	left_exons = unique(left_exons);
	right_exons = unique(right_exons);
	
	for ex = left_exons(1:end-1), fprintf(fid, '%s, ', ex{1}); end
	fprintf(fid, '%s\t', left_exons{end});
	
	for ex = right_exons(1:end-1), fprintf(fid, '%s, ', ex{1}); end
	fprintf(fid, '%s\t', right_exons{end});
	
	
	for k = 1:length(junction_types)-1
		fprintf(fid, '%s, ', junction_types{k});
	end
	fprintf(fid, '%s\t', junction_types{end});
	
	if ~isempty(frameshifts)
		for k = 1:length(frameshifts)-1
			fprintf(fid, '%d, ', frameshifts(k));
		end
		fprintf(fid, '%d\t', frameshifts(end));
	else
		fprintf(fid, '\t');
	end
	
	
	
	
	
	
	total_samples = [];
	for p = 1:num_exon_pairs
		total_samples = cat(1, total_samples, gfusion.ReadSamples{p});
	end
	
	% Display the total number of fusion positive patients at the end of the
	% spreadsheet.
	fprintf(fid, '%d\t', length(unique(total_samples)));
	
	if ~isempty(groups)
		for k = 1:2:length(groups)
			if ~isnumeric(groups{k+1}), error 'Invalid group format.'; end
			fprintf(fid, '%.1f\t', ...
				sum(ismember(total_samples, groups{k+1}))/length(groups{k+1}));
		end
	end

	
	% If we had associated gene or exon expression data available for the
	% samples, then show the logratios and p-values.
	if ~isempty(gene_expr)
		fprintf(fid, '%.2f\t%.2f\t', diff_expr(idx, 1), diff_expr(idx, 2));
	elseif ~isempty(exon_expr)
		fprintf(fid, '%.2f\t%.2f\t%.2f\t%.2f\t', diff_expr(idx, 1), ...
			diff_expr(idx, 2), diff_expr(idx, 3), diff_expr(idx, 4));
	end
	
	for s = 1:S
		sreads = sum(total_samples == s);
		if sreads == 0
			fprintf(fid, '\t');
		else
			fprintf(fid, '%d\t', sreads);
		end
	end
	
	fprintf(fid, '\n');
end

fclose(fid);

