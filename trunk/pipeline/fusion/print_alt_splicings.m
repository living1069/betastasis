function [] = print_alt_splicings(splices, varargin)

global organism;
genes = organism.Genes;
exons = organism.Exons;

min_reads = 0;
min_anchor_len = 10;
blacklist = {};
only_unique = false;
report_exon_seq = false;

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
	
	if strcmpi(varargin{k}, 'OnlyUnique')
		only_unique = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'ReportExonSeq')
		report_exon_seq = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isstruct(splices), splices = { splices }; end

S = length(splices);

joint_splices = {};
total_splice_reads = [];

% First we build a map of all gene pairs shown to participate in splices.
splice_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
for s = 1:S
	for f = 1:length(splices{s}.ReadCount)
		left_exon = splices{s}.Exons(f, 1);
		right_exon = splices{s}.Exons(f, 2);
		
		gene = organism.Exons.Gene(left_exon);
		
		key = gene;
		if ~splice_map.isKey(key)
			gsplice.Gene = gene;
			gsplice.Exons = [left_exon, right_exon];
			gsplice.ReadSequences = ...
				{ splices{s}.ReadSequences(f, 1:splices{s}.ReadCount(f))' };
			gsplice.JunctionOffsets = ...
				{ splices{s}.JunctionOffsets(f, 1:splices{s}.ReadCount(f)) };
			gsplice.ReadSamples = { repmat(s, splices{s}.ReadCount(f), 1) };
			joint_splices{end+1, 1} = gsplice;
			total_splice_reads(end+1, 1) = splices{s}.ReadCount(f);
			splice_map(key) = length(joint_splices);
		else
			fus_idx = splice_map(key);
			gsplice = joint_splices{fus_idx};
			idx = find(gsplice.Exons(:, 1) == left_exon & ...
				gsplice.Exons(:, 2) == right_exon);
			if isempty(idx)
				gsplice.Exons(end+1, :) = [left_exon, right_exon];
				gsplice.ReadSequences{end+1, 1} = ...
					splices{s}.ReadSequences(f, 1:splices{s}.ReadCount(f))';
				gsplice.JunctionOffsets{end+1, 1} = ...
					splices{s}.JunctionOffsets(f, 1:splices{s}.ReadCount(f));
				gsplice.ReadSamples{end+1, 1} = ...
					repmat(s, splices{s}.ReadCount(f), 1);
			else
				gsplice.ReadSequences{idx} = cat(1, ...
					gsplice.ReadSequences{idx}, ...
					splices{s}.ReadSequences(f, 1:splices{s}.ReadCount(f))');
				gsplice.JunctionOffsets{idx} = [ ...
					gsplice.JunctionOffsets{idx}, ...
					splices{s}.JunctionOffsets(f, 1:splices{s}.ReadCount(f)) ];
				gsplice.ReadSamples{idx} = [ gsplice.ReadSamples{idx}; ...
					repmat(s, splices{s}.ReadCount(f), 1) ];
			end
			total_splice_reads(fus_idx) = total_splice_reads(fus_idx) + ...
				splices{s}.ReadCount(f);
			joint_splices{fus_idx} = gsplice;
		end
	end
end

[~, order] = sort(total_splice_reads, 'descend');

fprintf(1, 'Potential splice transcripts (ordered by prevalence):\n');
for k = 1:length(order)
	idx = order(k);
	
	gsplice = joint_splices{idx};

	if total_splice_reads(idx) < min_reads, continue, end
	
	blacklisted = false;
	for k = 1:length(blacklist)
		if regexp(genes.Name{gsplice.Gene}, blacklist{k})
			blacklisted = true; break;
		end
	end
	
	if blacklisted, continue, end
		
	% Check the minimum anchor length.
	num_exon_pairs = size(gsplice.Exons, 1);
	best_anchor_len = 0;
	for p = 1:num_exon_pairs
		for s = 1:length(gsplice.ReadSequences{p})
			left_len = gsplice.JunctionOffsets{p}(s);
			seq = gsplice.ReadSequences{p}{s};
			right_len = length(seq) - left_len;
			
			best_anchor_len = max(min(left_len, right_len), best_anchor_len);
		end
	end

	if best_anchor_len < min_anchor_len, continue, end
	
	fprintf(1, '- splice variant of %s (%d total reads):\n', ...
		genes.Name{gsplice.Gene}, total_splice_reads(idx));
	
	num_exon_pairs = size(gsplice.Exons, 1);
	for p = 1:num_exon_pairs
		left_exon = gsplice.Exons(p, 1);
		right_exon = gsplice.Exons(p, 2);
		left_exon_seq = upper(exons.Sequence{left_exon});
		right_exon_seq = upper(exons.Sequence{right_exon});

		fprintf(1, '  * junction between %s[%s] and %s[%s]:\n', ...
			genes.Name{gsplice.Gene}, exons.ID{left_exon}, ...
			genes.Name{gsplice.Gene}, exons.ID{right_exon});
		
		for s = 1:length(gsplice.ReadSequences{p})
			left_len = gsplice.JunctionOffsets{p}(s);
			seq = gsplice.ReadSequences{p}{s};
			
			ref_seq = [left_exon_seq(end-left_len+1:end) ...
				right_exon_seq(1:length(seq)-left_len)];
			
			mismatches = (seq ~= ref_seq);
			seq(mismatches) = lower(seq(mismatches));

			fprintf(1, '    * %d: %s|', gsplice.ReadSamples{p}(s), ...
				seq(1:left_len));
			fprintf(1, '%s\n', seq(left_len+1:end));
		end
	end
end




