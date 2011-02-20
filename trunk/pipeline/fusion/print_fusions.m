function [] = print_fusions(fusions, varargin)

global organism;
genes = organism.Genes;
exons = organism.Exons;
chromosomes = organism.Chromosomes;

min_reads = 0;
min_anchor_len = 10;
blacklist = {};
only_unique = false;
report_exon_seq = false;

% Organism specific genes that are blacklisted by default.
if strcmpi(organism.Name, 'Homo sapiens')
	blacklist = { ...
		'RPPH1', 'LOC100008588', 'LOC100008589', 'RN18S1', 'SNORD', 'SNORA', ...
		'RNY', 'RN7SL', 'RN7SK', 'RNU', 'HLA-', ...
	};
end

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MinReads')
		min_reads = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Blacklist')
		blacklist = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MinAnchorLength')
		min_anchor_len = varargin{k+1};
		continue;
	end
	
	%if strcmpi(varargin{k}, 'OnlyUnique')
	%	only_unique = varargin{k+1};
	%	continue;
	%end
	
	if strcmpi(varargin{k}, 'ReportExonSeq')
		report_exon_seq = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isstruct(fusions)
	fusions = { fusions };
end

S = length(fusions);

joint_fusions = {};
total_fusion_reads = [];

print_sequences = false;
if isfield(fusions{1}, 'ReadSequences')
	print_sequences = true;
end

% First we build a map of all gene pairs shown to participate in fusions.
fusion_map = containers.Map;
for s = 1:S
	for f = 1:length(fusions{s}.ReadCount)
		left_exon = fusions{s}.Exons(f, 1);
		right_exon = fusions{s}.Exons(f, 2);
		
		left_gene = organism.Exons.Gene(fusions{s}.Exons(f, 1));
		right_gene = organism.Exons.Gene(fusions{s}.Exons(f, 2));
		
		key = sprintf('%d,%d', left_gene, right_gene);
		if ~fusion_map.isKey(key)
			gfusion.Genes = [left_gene, right_gene];
			gfusion.Exons = [left_exon, right_exon];
			gfusion.ReadSamples = { repmat(s, fusions{s}.ReadCount(f), 1) };
			
			if print_sequences
				gfusion.ReadSequences = ...
					{ fusions{s}.ReadSequences(f, 1:fusions{s}.ReadCount(f))' };
				gfusion.JunctionOffsets = ...
					{ fusions{s}.JunctionOffsets(f, 1:fusions{s}.ReadCount(f))};
			end
			
			joint_fusions{end+1, 1} = gfusion;
			total_fusion_reads(end+1, 1) = fusions{s}.ReadCount(f);
			fusion_map(key) = length(joint_fusions);
		else
			fus_idx = fusion_map(key);
			gfusion = joint_fusions{fus_idx};
			R = fusions{s}.ReadCount(f);
			
			idx = find(gfusion.Exons(:, 1) == left_exon & ...
				gfusion.Exons(:, 2) == right_exon);
			if isempty(idx)
				gfusion.Exons(end+1, :) = [left_exon, right_exon];
				gfusion.ReadSamples{end+1, 1} = repmat(s, R, 1);

				if print_sequences
					gfusion.ReadSequences{end+1, 1} = ...
						fusions{s}.ReadSequences(f, 1:R)';
					gfusion.JunctionOffsets{end+1, 1} = ...
						fusions{s}.JunctionOffsets(f, 1:R);
				end
			else
				gfusion.ReadSamples{idx} = [ gfusion.ReadSamples{idx}; ...
					repmat(s, R, 1) ];

				if print_sequences
					gfusion.ReadSequences{idx} = cat(1, ...
						gfusion.ReadSequences{idx}, ...
						fusions{s}.ReadSequences(f, 1:R)');
					gfusion.JunctionOffsets{idx} = [ ...
						gfusion.JunctionOffsets{idx}, ...
						fusions{s}.JunctionOffsets(f, 1:R) ];
				end
			end
			total_fusion_reads(fus_idx) = total_fusion_reads(fus_idx) + R;
			joint_fusions{fus_idx} = gfusion;
		end
	end
end

if only_unique
	for j = 1:length(joint_fusions)
		gfusion = joint_fusions{j};
		%unique_seq = unique(gfusion.ReadSequences
	end
end





% Discard fusions according to a number of possible criteria.
keep = false(size(total_fusion_reads));
for f = 1:length(total_fusion_reads)
	gfusion = joint_fusions{f};

	if total_fusion_reads(f) < min_reads, continue, end
	
	blacklisted = false;
	for k = 1:length(blacklist)
		if regexp(genes.Name{gfusion.Genes(1)}, blacklist{k})
			blacklisted = true; break;
		end
		if regexp(genes.Name{gfusion.Genes(2)}, blacklist{k})
			blacklisted = true; break;
		end
	end
	
	if blacklisted, continue, end
		
	% Check the minimum anchor length.
	if isfield(gfusion, 'ReadSequences')
		num_exon_pairs = size(gfusion.Exons, 1);
		best_anchor_len = 0;
		for p = 1:num_exon_pairs
			for s = 1:length(gfusion.ReadSequences{p})
				left_len = gfusion.JunctionOffsets{p}(s);
				seq = gfusion.ReadSequences{p}{s};
				right_len = length(seq) - left_len;
				
				best_anchor_len = max( ...
					min(left_len, right_len), best_anchor_len);
			end
		end
		
		if best_anchor_len < min_anchor_len, continue, end
	end
	
	keep(f) = true;    % Fusion passed all criteria.
end

fprintf(1, '%d / %d (%.1f%%) fusions passed all criteria.\n', sum(keep), ...
	length(total_fusion_reads), sum(keep) / length(total_fusion_reads) * 100);
total_fusion_reads = total_fusion_reads(keep);
joint_fusions = joint_fusions(keep);





% Print all fusions, starting with the ones with the highest level of evidence.
if 0
	[~, order] = sort(total_fusion_reads, 'descend');
else
	chisq = nan(size(total_fusion_reads));
	for f = 1:length(total_fusion_reads)
		gfusion = joint_fusions{f};
		num_exon_pairs = size(gfusion.Exons, 1);
		
		sr = zeros(1, S);
		for p = 1:num_exon_pairs
			for s = 1:S, sr(s) = sum(gfusion.ReadSamples{p} == s); end
		end
		expected_reads = sum(sr) / S;
		chisq(f) = sum((sr - ones(1, S) * expected_reads).^2 / expected_reads);
	end
	
	[~, order] = sort(chisq, 'descend');
end

fprintf(1, 'Potential fusion transcripts (ordered by prevalence):\n');
for k = 1:length(order)
	idx = order(k);
	gfusion = joint_fusions{idx};
	
	fprintf(1, '- fusion of %s and %s (%d total reads):\n', ...
		genes.Name{gfusion.Genes(1)}, genes.Name{gfusion.Genes(2)}, ...
		total_fusion_reads(idx));
	
	num_exon_pairs = size(gfusion.Exons, 1);
	for p = 1:num_exon_pairs
		left_exon = gfusion.Exons(p, 1);
		right_exon = gfusion.Exons(p, 2);
		left_exon_seq = upper(exons.Sequence{left_exon});
		right_exon_seq = upper(exons.Sequence{right_exon});
		
		fprintf(1, '  * junction between %s[%s] and %s[%s]:\n', ...
			genes.Name{gfusion.Genes(1)}, exons.ID{gfusion.Exons(p, 1)}, ...
			genes.Name{gfusion.Genes(2)}, exons.ID{gfusion.Exons(p, 2)});
		
		if report_exon_seq
			left_chr = genes.Chromosome(exons.Gene(left_exon));
			right_chr = genes.Chromosome(exons.Gene(right_exon));
			
			if isnan(left_chr)
				left_chr = 'unknown chr';
			else
				left_chr = ['chr' chromosomes.Name{left_chr}];
			end
			
			if isnan(right_chr)
				right_chr = 'unknown chr';
			else
				right_chr = ['chr' chromosomes.Name{right_chr}];
			end
			
			left_strand = genes.Strand(exons.Gene(left_exon));
			right_strand = genes.Strand(exons.Gene(right_exon));
			
			if left_strand == ' ', left_strand = 'unknown'; end
			if right_strand == ' ', right_strand = 'unknown'; end
			
			fprintf(1, '    * %s[%s] sequence (%s, %s strand): %s\n', ...
				genes.Name{exons.Gene(gfusion.Exons(p, 1))}, ...
				exons.ID{gfusion.Exons(p, 1)}, left_chr, left_strand, ...
				left_exon_seq);
			fprintf(1, '    * %s[%s] sequence (%s, %s strand): %s\n', ...
				genes.Name{exons.Gene(gfusion.Exons(p, 2))}, ...
				exons.ID{gfusion.Exons(p, 2)}, right_chr, right_strand, ...
				right_exon_seq);
		end
		
		if print_sequences
			for s = 1:length(gfusion.ReadSequences{p})
				left_len = gfusion.JunctionOffsets{p}(s);
				seq = gfusion.ReadSequences{p}{s};
				
				ref_seq = [left_exon_seq(end-left_len+1:end) ...
					right_exon_seq(1:length(seq)-left_len)];
				
				mismatches = (seq ~= ref_seq);
				seq(mismatches) = lower(seq(mismatches));
				
				fprintf(1, '    * %d: %s|', gfusion.ReadSamples{p}(s), ...
					seq(1:left_len));
				fprintf(1, '%s\n', seq(left_len+1:end));
			end
		else
			num_samples = max(gfusion.ReadSamples{p});
			for s = 1:num_samples
				sr = sum(gfusion.ReadSamples{p} == s);
				if sr == 0, continue, end
				fprintf(1, '    * reads in sample %d: %d\n', s, sr);
			end
		end
	end
end

