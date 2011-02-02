function [] = print_fusions(fusions, varargin)

global organism;
genes = organism.Genes;
exons = organism.Exons;

min_reads = 0;
min_anchor_len = 10;
blacklist = {};

% Organism specific genes that are blacklisted by default.
if strcmpi(organism.Name, 'Homo sapiens')
	blacklist = { ...
		'RPPH1', 'LOC100008588', 'LOC100008589', 'RN18S1', 'SNORD', 'SNORA', ...
		'RNY1', 'RNY3', 'RNY5', 'RN7SL', 'RN7SK', 'RNU2-1' 'RNU5D', ...
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
	
	error('Unrecognized option "%s".', varargin{k});
end

if isstruct(fusions)
	fusions = { fusions };
end

S = length(fusions);

joint_fusions = {};
total_fusion_reads = [];

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
			gfusion.ReadSequences = ...
				{ fusions{s}.ReadSequences(f, 1:fusions{s}.ReadCount(f))' };
			gfusion.JunctionOffsets = ...
				{ fusions{s}.JunctionOffsets(f, 1:fusions{s}.ReadCount(f)) };
			gfusion.ReadSamples = { repmat(s, fusions{s}.ReadCount(f), 1) };
			joint_fusions{end+1, 1} = gfusion;
			total_fusion_reads(end+1, 1) = fusions{s}.ReadCount(f);
			fusion_map(key) = length(joint_fusions);
		else
			fus_idx = fusion_map(key);
			gfusion = joint_fusions{fus_idx};
			idx = find(gfusion.Exons(:, 1) == left_exon & ...
				gfusion.Exons(:, 2) == right_exon);
			if isempty(idx)
				gfusion.Exons(end+1, :) = [left_exon, right_exon];
				gfusion.ReadSequences{end+1, 1} = ...
					fusions{s}.ReadSequences(f, 1:fusions{s}.ReadCount(f))';
				gfusion.JunctionOffsets{end+1, 1} = ...
					fusions{s}.JunctionOffsets(f, 1:fusions{s}.ReadCount(f));
				gfusion.ReadSamples{end+1, 1} = ...
					repmat(s, fusions{s}.ReadCount(f), 1);
			else
				gfusion.ReadSequences{idx} = cat(1, ...
					gfusion.ReadSequences{idx}, ...
					fusions{s}.ReadSequences(f, 1:fusions{s}.ReadCount(f))');
				gfusion.JunctionOffsets{idx} = [ ...
					gfusion.JunctionOffsets{idx}, ...
					fusions{s}.JunctionOffsets(f, 1:fusions{s}.ReadCount(f)) ];
				gfusion.ReadSamples{idx} = [ gfusion.ReadSamples{idx}; ...
					repmat(s, fusions{s}.ReadCount(f), 1) ];
			end
			total_fusion_reads(fus_idx) = total_fusion_reads(fus_idx) + ...
				fusions{s}.ReadCount(f);
			joint_fusions{fus_idx} = gfusion;
		end
	end
end

[~, order] = sort(total_fusion_reads, 'descend');

fprintf(1, 'Potential fusion transcripts (ordered by prevalence):\n');
for k = 1:length(order)
	idx = order(k);
	
	gfusion = joint_fusions{idx};

	if total_fusion_reads(idx) < min_reads, continue, end
	
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
	num_exon_pairs = size(gfusion.Exons, 1);
	best_anchor_len = 0;
	for p = 1:num_exon_pairs
		for s = 1:length(gfusion.ReadSequences{p})
			left_len = gfusion.JunctionOffsets{p}(s);
			seq = gfusion.ReadSequences{p}{s};
			right_len = length(seq) - left_len;
			
			best_anchor_len = max(min(left_len, right_len), best_anchor_len);
		end
	end
	
	if best_anchor_len < min_anchor_len, continue, end
		
	fprintf(1, '- fusion of %s and %s (%d total reads):\n', ...
		genes.Name{gfusion.Genes(1)}, genes.Name{gfusion.Genes(2)}, ...
		total_fusion_reads(idx));
	
	num_exon_pairs = size(gfusion.Exons, 1);
	for p = 1:num_exon_pairs
		fprintf(1, '  * junction between %s[%s] and %s[%s]:\n', ...
			genes.Name{gfusion.Genes(1)}, exons.ID{gfusion.Exons(p, 1)}, ...
			genes.Name{gfusion.Genes(2)}, exons.ID{gfusion.Exons(p, 2)});
			
		left_exon_seq = upper(exons.Sequence{gfusion.Exons(p, 1)});
		right_exon_seq = upper(exons.Sequence{gfusion.Exons(p, 2)});
		
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
	end
end

