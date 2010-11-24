function [] = print_snps(snps, threshold)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;

if nargin < 2, threshold = 0.1; end
	
[~, order] = sort(snps.Score, 1, 'descend');
fprintf(1, 'SNPs ordered by significance score:\n');
for k = 1:length(order)
	idx = order(k);
	if snps.Score(idx) < threshold, continue, end
	
	signature = snps.Signature{idx};
	pos = find(signature == ':');
	ts = signature(1:pos(1)-1);
	if regexpi(ts, '^(NM_|XM_|NR_|XR_).+')
		signature = [genes.Name{transcripts.Gene(transcript_idx(ts))} ...
			', ' signature];
	else
		% Filter out microRNA mutations where the nucleotide change occurs
		% at the 3' end.
		if str2num(signature(pos(1)+1:pos(2)-1)) >= 15
			continue;
		end
	end
	
	
	fprintf(1, '- %s: %d/%d (%.1f%%)\n', signature, ...
		snps.SupportingReads(idx), snps.TotalReads(idx), ...
		snps.SupportingReads(idx) / snps.TotalReads(idx) * 100);
end

