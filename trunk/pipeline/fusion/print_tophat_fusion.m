
% Author: Matti Annala <matti.annala@tut.fi>

function [] = print_tophat_fusions()

global organism;
genes = organism.Genes;

min_distance = 1e6;
min_longest_anchor = 30;
min_single_reads = 5;
min_paired_reads = 0;




% Now we parse the fusions.out results and load them into Matlab.
fprintf(1, 'Parsing fusion candidates...\n');

files = find_files('fusions.out');
samples = regexprep(files, '.*/([^/]+)/tophat_out/fusions.out', '$1');

S = length(samples);

% These variables are used for collecting fusion statistics across all samples.
fusion_names = {};
fusion_reads = zeros(0, S);
fusion_map = containers.Map;

for s = 1:S
	fprintf(1, 'Reporting fusion candidates for sample %s...\n', samples{s});
	
	fid = fopen(files{s});
	data = textscan(fid, ...
		'%s%d%d%s%d%d%d%d%d%d %*s%*s%*s%*s %s %*s %s %*s %s %*s %s %*s%*s', ...
		'Delimiter', '\t', 'ReturnOnError', false);
	fclose(fid);
	
	chrs = data{1};
	tokens = regexp(chrs, '(.+)-(.+)', 'tokens');
	
	for k = 1:length(chrs)
		left_chr{k} = tokens{k}{1}{1};
		right_chr{k} = tokens{k}{1}{2};
	end
	
	left_chr_num = chromosome_sym2num(left_chr);
	right_chr_num = chromosome_sym2num(right_chr);
	
	left_offset = data{2};
	right_offset = data{3};
	orientation = data{4};
	single_span_reads = data{5};
	paired_span_reads = data{6};
	combo_span_reads = data{7};
	%contradicting_reads = data{8};
	longest_left_anchor = data{9};
	longest_right_anchor = data{10};
	left_contig = data{11};
	right_contig = data{12};
	%left_evidence = data{13};
	%right_evidence = data{14};
	
	
	% Score solely based on junction spanning single end reads.
	%score = min(single_span_reads, paired_span_reads);
	score = single_span_reads;
	[~, order] = sort(score, 'descend');
	
	fid = fopen(sprintf('fusion_candidates_%s.txt', samples{s}), 'W');
	fprintf(fid, ['LEFT_CHR\tRIGHT_CHR\tLEFT_OFFSET\tRIGHT_OFFSET\t' ...
		'LEFT_GENES\tRIGHT_GENES\tLEFT_SEQ\tRIGHT_SEQ\t' ...
		'SINGLE + PAIRED + COMBO\n']);
		
	for k = order'
		if longest_left_anchor(k) < min_longest_anchor || ...
			longest_right_anchor(k) < min_longest_anchor || ...
			single_span_reads(k) < min_single_reads || ...
			paired_span_reads(k) < min_paired_reads
			continue;
		end
		
		if left_chr_num(k) == right_chr_num(k) && ...
			abs(left_offset(k) - right_offset(k)) < min_distance
			continue;
		end
		
		% Check if the 5' sides or 3' sides of the breakpoints exhibit a large
		% degree of similarity. If they do, discard the fusion as a likely
		% sequencing artifact.
		al_score = swalign(left_contig{k}(1:50), right_contig{k}(1:50));
		if al_score > 90, continue, end

		al_score = swalign(left_contig{k}(52:101), right_contig{k}(52:101));
		if al_score > 90
			%fprintf(1, 'Discarded %s\n          %s\n', ...
			%	left_contig{k}(52:101), right_contig{k}(52:101));
			continue;
		end
		
		% FIXME: Check if the junction sequence as a whole aligns to some
		% contiguous location in the genome. If this happens, then we have
		% to assume that the read actually originates from this location.
		
		
		fprintf(fid, '%s\t%s\t%d\t%d\t', ...
			left_chr{k}, right_chr{k}, left_offset(k), right_offset(k));
			
		% Calculate which genes overlap with the fusion breakpoints on the left
		% and right side.
		left_overlap_genes = find(genes.Chromosome == left_chr_num(k) & ...
			genes.Position(:, 1) - 10e3 <= left_offset(k) & ...
			genes.Position(:, 2) + 10e3 >= left_offset(k))';
		right_overlap_genes = find(genes.Chromosome == right_chr_num(k) & ...
			genes.Position(:, 1) - 10e3 <= right_offset(k) & ...
			genes.Position(:, 2) + 10e3 >= right_offset(k))';
			
		for l = left_overlap_genes'
			for r = right_overlap_genes'
				key = sprintf('%s-%s', genes.Name{l}, genes.Name{r});
				if ~fusion_map.isKey(key)
					fusion_reads(end+1, 1) = 0;
					fusion_map(key) = size(fusion_reads, 1);
					fusion_names{size(fusion_reads, 1)} = key;
				end
				fusion_reads(fusion_map(key), s) = ...
					fusion_reads(fusion_map(key), s) + single_span_reads(k);
			end
		end
		
		fprintf_cell(fid, genes.Name(left_overlap_genes), ', ');
		fprintf(fid, '\t');
		fprintf_cell(fid, genes.Name(right_overlap_genes), ', ');
		
		fprintf(fid, '\t%s|%s\t%d + %d + %d\n', ...
			left_contig{k}(1:50), right_contig{k}(52:101), ...
			single_span_reads(k), paired_span_reads(k), combo_span_reads(k));
	end
	fclose(fid);
end


% Print a list of the most promising fusions based on chi-square test.
chisq = nan(1, size(fusion_reads, 1));
for f = 1:length(chisq)
	expected_reads = sum(fusion_reads(f, :)) / S;
	chisq(f) = sum((fusion_reads(f, :) - ...
		ones(1, S) * expected_reads).^2 / expected_reads);
end
	
[~, order] = sort(chisq, 'descend');

fid = fopen('best_fusions_chisq.txt', 'W');
fprintf(fid, 'FUSION');
fprintf(fid, '\t%s', samples{:});
%for s = 1:S, fprintf(fid, '\t%s', samples{s}); end
fprintf(fid, '\n');
for k = order
	fprintf(fid, '%s', fusion_names{k});
	fprintf(fid, '\t%d', fusion_reads(k, :));
	fprintf(fid, '\n');
end
fclose(fid);





