function joint_fusions = filter_fusions(rearrangements, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;
chromosomes = organism.Chromosomes;

min_frequency = 0;
min_distance = 0;
min_anchor_len = 0;
blacklist = {};
min_read_anchor_len = 0;
negative_samples = [];
max_avg_mismatches = Inf;
min_homology_distance = 0;
whitelist = {};

% Organism specific genes that are blacklisted by default.
if strcmpi(organism.Name, 'Homo sapiens')
	blacklist = { ...
		'RPPH1', 'LOC100008588', 'LOC100008589', 'RN18S1', 'SNORD', 'SNORA', ...
		'RNY', 'RN7SL', 'RN7SK', 'RNU', 'HLA-', ...
	};
end

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'whitelist')
		whitelist = varargin{k+1};
		if ischar(whitelist), whitelist = { whitelist }; end
		continue;
	end
	
	if rx(varargin{k}, 'blacklist')
		blacklist = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'min.*freq')
		min_frequency = varargin{k+1}; continue;
	end
	
	if regexpi(varargin{k}, 'max.*(avg|average).*mismatch')
		max_avg_mismatches = varargin{k+1}; continue;
	end
	
	if strcmpi(varargin{k}, 'MinDistance')
		min_distance = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'MinHomologyDistance')
		min_homology_distance = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'MinAnchorLen')
		min_anchor_len = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'MinReadAnchorLen')
		min_read_anchor_len = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'NegativeSamples')
		negative_samples = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

fusions = rearrangements.fusions;
S = length(fusions);

joint_fusions = {};


% First we build a map of all gene pairs shown to participate in fusions.
fusion_map = containers.Map;
for s = 1:S

	for f = 1:length(fusions{s}.ReadCount)
		left_exon = fusions{s}.Exons(f, 1);
		right_exon = fusions{s}.Exons(f, 2);
		
		left_gene = exons.Gene(fusions{s}.Exons(f, 1));
		right_gene = exons.Gene(fusions{s}.Exons(f, 2));
		
		key = sprintf('%d,%d', left_gene, right_gene);
		if ~fusion_map.isKey(key)
			gfusion.Genes = [left_gene, right_gene];
			gfusion.Exons = [left_exon, right_exon];
			gfusion.ReadSamples = { repmat(s, fusions{s}.ReadCount(f), 1) };
			
			gfusion.ReadSequences = { fusions{s}.ReadSequences{f} };
			gfusion.JunctionOffsets = { fusions{s}.JunctionOffsets{f} };
			
			joint_fusions{end+1, 1} = gfusion;
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

				gfusion.ReadSequences{end+1, 1} = fusions{s}.ReadSequences{f};
				gfusion.JunctionOffsets{end+1, 1} = ...
					fusions{s}.JunctionOffsets{f};
			else
				gfusion.ReadSamples{idx} = [ gfusion.ReadSamples{idx}; ...
					repmat(s, R, 1) ];

				gfusion.ReadSequences{idx} = cat(1, ...
					gfusion.ReadSequences{idx}, fusions{s}.ReadSequences{f});
				gfusion.JunctionOffsets{idx} = [ ...
					gfusion.JunctionOffsets{idx}; ...
					fusions{s}.JunctionOffsets{f} ];
			end
			joint_fusions{fus_idx} = gfusion;
		end
	end
end


total_fusion_reads = 0;
total_fusion_juncs = 0;
for f = 1:length(joint_fusions)
	total_fusion_juncs = total_fusion_juncs + size(joint_fusions{f}.Exons, 1);
	total_fusion_reads = total_fusion_reads + ...
		sum(cellfun(@length, joint_fusions{f}.ReadSamples));
end

fprintf('%d total reads aligned to %d fusion junctions...\n', ...
	total_fusion_reads, total_fusion_juncs);




fprintf('Filtering fusion candidates according to specified criteria...\n');
progress = Progress;

% Discard fusions according to a number of criteria.
discard = zeros(size(joint_fusions));
for f = 1:length(joint_fusions)
	progress.update(f / length(joint_fusions));
	gfusion = joint_fusions{f};
	
	% If the user defined a whitelist, check that fusion gene name matches
	% with the regexp.
	if ~isempty(whitelist)
		whitelisted = false;
		for k = 1:length(whitelist)
			name = [genes.Name{gfusion.Genes(1)} '-' ...
				genes.Name{gfusion.Genes(2)}];
			if regexp(name, whitelist{k})
				whitelisted = true; break;
			end
		end
		if ~whitelisted, discard(f) = 1; continue; end
	end

	% Discard fusions based on a blacklist of genes that seem to be "fused"
	% all over the place. The default blacklist includes rRNA genes and similar
	% very highly expressed genes.
	blacklisted = false;
	for k = 1:length(blacklist)
		if regexp(genes.Name{gfusion.Genes(1)}, blacklist{k})
			blacklisted = true; break;
		end
		if regexp(genes.Name{gfusion.Genes(2)}, blacklist{k})
			blacklisted = true; break;
		end
	end
	if blacklisted, discard(f) = 1; continue; end
		
	% Check the minimum anchor length.
	if isfield(gfusion, 'ReadSequences')
		num_exon_pairs = size(gfusion.Exons, 1);
		best_anchor_len = [0 0];
		
		for p = 1:num_exon_pairs
			% We also discard all reads that do not exceed the minimum
			% per-read anchor length.
			keep = true(1, length(gfusion.ReadSequences{p}));
			
			for s = 1:length(gfusion.ReadSequences{p})
				left_len = gfusion.JunctionOffsets{p}(s);
				right_len = length(gfusion.ReadSequences{p}{s}) - left_len;
				
				if min(left_len, right_len) < min_read_anchor_len
					keep(s) = false;
				end
				
				best_anchor_len = max(best_anchor_len, [left_len right_len]);
			end
			
			if min_read_anchor_len > 0
				gfusion.ReadSamples{p} = gfusion.ReadSamples{p}(keep);
				gfusion.ReadSequences{p} = gfusion.ReadSequences{p}(keep);
				gfusion.JunctionOffsets{p} = gfusion.JunctionOffsets{p}(keep);
			end
		end
		
		if min_read_anchor_len > 0
			keep = ~cellfun(@isempty, gfusion.ReadSamples);
			if all(~keep)
				discard(f) = 5; continue;
			else
				gfusion.Exons = gfusion.Exons(keep, :);
				gfusion.ReadSamples = gfusion.ReadSamples(keep);
				gfusion.ReadSequences = gfusion.ReadSequences(keep);
				gfusion.JunctionOffsets = gfusion.JunctionOffsets(keep);
				joint_fusions{f} = gfusion;
			end
		end
		
		if min(best_anchor_len) < min_anchor_len
			discard(f) = 5; continue;
		end
	end
		
	if min_frequency > 0
		% Make it possible to specify frequency as number of samples.
		if min_frequency < 1, min_frequency = min_frequency * S; end
		
		samples_with_fusion = false(S, 1);
		for p = 1:size(gfusion.Exons, 1)
			samples_with_fusion(gfusion.ReadSamples{p}) = true;
		end
		
		if sum(samples_with_fusion) < min_frequency
			discard(f) = 2; continue;
		end
	end
	
	% Discard fusions that show reads in known fusion negative samples.
	if ~isempty(negative_samples)
		samples_with_fusion = false(S, 1);
		for p = 1:size(gfusion.Exons, 1)
			samples_with_fusion(gfusion.ReadSamples{p}) = true;
		end
		
		if any(samples_with_fusion(negative_samples))
			discard(f) = 3; continue;
		end
	end
	
	if min_distance > 0
		chr = genes.Chromosome(gfusion.Genes);
		pos = genes.Position(gfusion.Genes, :);
		
		if chr(1) == chr(2) && ...
			abs(mean(pos(1, :)) - mean(pos(2, :))) < min_distance
			discard(f) = 4; continue;
		end
	end
		
	% Check that the reads do not have an aberrantly high number of mismatches.
	if max_avg_mismatches < Inf
		if isfield(gfusion, 'ReadSequences')
			mismatches = [0 0];
			for p = 1:size(gfusion.Exons, 1)
				for k = 1:length(gfusion.ReadSequences{p})
					left_exon_seq = upper(exons.Sequence{gfusion.Exons(p, 1)});
					right_exon_seq = upper(exons.Sequence{gfusion.Exons(p, 2)});
					
					left_len = gfusion.JunctionOffsets{p}(k);
					seq = gfusion.ReadSequences{p}{k};
					
					ref_seq = [left_exon_seq(end-left_len+1:end) ...
						right_exon_seq(1:length(seq)-left_len)];
					
					mismatches = mismatches + [sum(seq ~= ref_seq), 1];
				end
			end
			
			if mismatches(1) / mismatches(2) > max_avg_mismatches
				discard(f) = 6; continue;
			end
		end
	end
	
	% Check if either side of the fusion junction is found in the genomic
	% vicinity of the other gene.
	if min_homology_distance > 0
		if ~fusion_check_neighborhood( ...
			[genes.Name{gfusion.Genes(1)} ':' genes.Name{gfusion.Genes(2)}], ...
			gfusion, min_homology_distance)
			discard(f) = 7; continue;
		end
	end
end

fprintf('%d / %d (%.1f%%) discarded due to blacklisted genes.\n', ...
	sum(discard == 1), length(discard), ...
	100 * sum(discard == 1) / length(discard));
fprintf('%d / %d (%.1f%%) discarded due to low frequency.\n', ...
	sum(discard == 2), length(discard), ...
	100 * sum(discard == 2) / length(discard));
fprintf('%d / %d (%.1f%%) discarded due to presence in controls.\n', ...
	sum(discard == 3), length(discard), ...
	100 * sum(discard == 3) / length(discard));
fprintf('%d / %d (%.1f%%) discarded due to too close proximity.\n', ...
	sum(discard == 4), length(discard), ...
	100 * sum(discard == 4) / length(discard));
fprintf('%d / %d (%.1f%%) discarded due to inadequate anchor length.\n', ...
	sum(discard == 5), length(discard), ...
	100 * sum(discard == 5) / length(discard));
fprintf('%d / %d (%.1f%%) discarded due to recurrent mismatches.\n', ...
	sum(discard == 6), length(discard), ...
	100 * sum(discard == 6) / length(discard));
fprintf('%d / %d (%.1f%%) discarded due to genomic neighborhood.\n', ...
	sum(discard == 7), length(discard), ...
	100 * sum(discard == 7) / length(discard));


fprintf('%d / %d (%.1f%%) fusions passed all criteria.\n', ...
	sum(discard == 0), length(joint_fusions), ...
	sum(discard == 0) / length(joint_fusions) * 100);
joint_fusions = joint_fusions(discard == 0);


