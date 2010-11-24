function pval = snp_significance(sample_mutations, ref_mutations)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;

signatures = {};
signature_to_idx = containers.Map;
ratios = [];
supporting_reads = [];

mutation_snps = cat(2, sample_mutations.SNPs, ref_mutations.SNPs);

num_sigs = 0;
for k = 1:length(mutation_snps)
	snps = mutation_snps{k};
	for s = 1:length(snps.Signature)
		sig = snps.Signature{s};
		if ~signature_to_idx.isKey(sig)
			num_sigs = num_sigs + 1;
			signature_to_idx(sig) = num_sigs;
			signatures{end + 1} = sig;
		end
		
		idx = signature_to_idx(sig);
		ratios(idx, k) = snps.SupportingReads(s) / snps.TotalReads(s);
		supporting_reads(idx, k) = snps.SupportingReads(s);
	end
end

sample_indices = 1:length(sample_mutations.SNPs);
ref_indices = length(sample_mutations.SNPs)+1:size(ratios, 2);

[~, p] = ttest2(ratios(:, sample_indices)', ratios(:, ref_indices)');
pval = p';

fprintf(1, 'SNPs with the most significant difference between categories:\n');
[~, order] = sort(pval);
for k = 1:length(order)
	idx = order(k);
	if pval(idx) >= 0.05 || isnan(pval(idx)), continue, end
		
	signature = signatures{idx};
	pos = find(signature == ':');
	ts = signature(1:pos(1)-1);
	if regexpi(ts, '^(NM_|XM_|NR_|XR_).+')
		signature = [genes.Name{transcripts.Gene(transcript_idx(ts))} ...
			', ' signature];
	end
		
	fprintf(1, '- %s: %.1f vs %.1f (p = %f, %d reads total)\n', ...
		signature, median(ratios(idx, sample_indices)) * 100, ...
		median(ratios(idx, ref_indices)) * 100, pval(idx), ...
		sum(supporting_reads(idx, :)));
end

