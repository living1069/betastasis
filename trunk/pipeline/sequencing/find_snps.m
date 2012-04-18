function mutations = find_snps(reads, transcripts, varargin)

global organism;

threshold = 0.05;
adapter = '';

S = length(reads.Raw);

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Adapter')
		adapter = upper(varargin{k+1});
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

adapter_len = length(adapter);

transcripts.Sequence = strcat(transcripts.Sequence, adapter);

transcript_name_to_idx = containers.Map(transcripts.Name, ...
	num2cell(1:length(transcripts.Name)));
	
mutations = struct;
mutations.SNPs = {};

mutations.Meta = reads.Meta;
mutations.Meta.Type = 'SNPs';
mutations.Meta.Organism = organism.Name;
mutations.Meta.ScoreThreshold = repmat([threshold], S, 1);

for s = 1:S
	fprintf(1, 'Analyzing sample %s for mutations:\n', ...
		reads.Meta.Sample.Filename{s});
	
	al = align_reads(filter_query(reads, s), transcripts, ...
		'MaxMismatches', 2, 'AllowAlignments', 20, ...
		'Columns', 'target,offset,sequence', varargin{:});
	
	fprintf(1, 'Calculating nucleotide frequencies...\n');
	tx_idx = cell2mat(transcript_name_to_idx.values(al.Target));
	
	% Sort the reads based on the transcripts to which they align.
	[~, order] = sort(tx_idx);
	al.Offset = al.Offset(order);
	al.Sequence = al.Sequence(order);
	
	% We preallocate the vectors to avoid triggering zillions of dynamic
	% memory allocations during the SNP discovery process.
	snps = struct('Signature', { cell(10000, 1) }, ...
	              'Score', zeros(10000, 1), ...
				  'SupportingReads', zeros(10000, 1), ...
				  'TotalReads', zeros(10000, 1));
	
	progress = Progress;
	found = 0;
	
	% In order to reduce memory usage, we partition the transcripts into sets
	% that are analyzed separately.
	tx = 1;
	nuc_reads = [];
	for k = 1:length(tx_idx)
		if tx_idx(k) ~= tx
			progress.update(k / length(tx_idx));
			
			if ~isempty(nuc_reads)
				totals = sum(nuc_reads, 1);
				seq = transcripts.Sequence{k};
				nuc_nums = nt2int(seq);
				
				scores = tanh(nuc_reads/100) .* nuc_reads ./ repmat(totals,4,1);
				
				for s = 1:length(seq)
					for b = setdiff(1:4, nuc_nums(s))
						%if nuc_reads(b, s) < 10, continue, end

						if scores(b, s) > threshold
							found = found + 1;
							snps.Signature{found} = sprintf('%s:%d:%s>%s', ...
								transcripts.Name{tx}, s, seq(s), int2nt(b));
							snps.Score(found) = scores(b, s);
							snps.SupportingReads(found) = nuc_reads(b, s);
							snps.TotalReads(found) = totals(s);
						end
					end
				end
			end
			
			tx = tx_idx(k);
			nuc_reads = zeros(4, length(transcripts.Sequence{tx}) + ...
				adapter_len + 20, 'single');
		end
		
		% Linear indexing is used to efficiently increment the matrix
		% elements  that correspond to the read sequence bases.
		pos = ((al.Offset(k): ...
			al.Offset(k) + length(al.Sequence{k}) - 1) - 1) * 4 + ...
			int32(nt2int(al.Sequence{k}));
		
		nuc_reads(pos) = nuc_reads(pos) + 1;
	end
		
		%fprintf(1, 'Calculating background probabilities of read errors...\n');
		%snp_probs = zeros(4, 4);
		%for orig = 1:4
		%	for mod = setdiff(1:4, orig)
		%		total_reads = 0;
		%		total_mod_nuc = 0;
		%		for k = 1:length(transcripts.Name)
		%			freqs = nuc_reads{k};
		%			if isempty(freqs), continue, end
		%			
		%			seq = transcripts.Sequence{k};
		%			orig_nuc_idx = (seq == int2nt(orig));
		%			
		%			total_mod_nuc = total_mod_nuc + sum(freqs(mod, orig_nuc_idx));
		%			total_reads = total_reads + sum(sum(freqs(:, orig_nuc_idx)));
		%		end
		%		
		%		snp_probs(orig, mod) = total_mod_nuc / total_reads; 
		%	end
		%end
		
		%snp_probs
		%snp_probs = snp_probs / mean(mean(snp_probs));
		
	if 0
		offset_errors = zeros(1, 100);
		offset_reads = zeros(1, 100);
		
		offset_snp_ratio_sums = zeros(1, 100);
		length_cum = zeros(1, 100);
		
		for k = 1:length(transcripts.Name)
			reads = nuc_reads{k};
			if isempty(reads), continue, end
			
			seq = transcripts.Sequence{k};
			for s = 1:length(seq)
				misbases = setdiff(1:4, nt2int(seq(s)));
				offset_errors(s) = offset_errors(s) + sum(reads(misbases, s));
				offset_reads(s) = offset_reads(s) + sum(reads(:, s));
				
				if offset_reads(s) > 0
					offset_snp_ratio_sums(s) = offset_snp_ratio_sums(s) + ...
						offset_errors(s) / offset_reads(s);
				end
				length_cum(s) = length_cum(s) + 1;
			end
		end
		
		offset_probs = offset_errors ./ offset_reads;
		offset_probs(isnan(offset_probs)) = 0;
	
		figure; bar(offset_probs);
		xlim([0 25.5]);
		xlabel('Nucleotide offset'); ylabel('Probability of read error');
		saveas(gcf, '~/mirna_snp_rate.pdf');
		
		foo = offset_snp_ratio_sums ./ length_cum;
		foo(isnan(foo)) = 0;
		
		figure; bar(foo);
		xlim([0 25.5]);
		xlabel('Nucleotide offset'); ylabel('Average of SNP rates');
		saveas(gcf, '~/mirna_snp_prob_avg.pdf');
		
		offset_probs = offset_errors ./ offset_reads;
		offset_probs = offset_probs / nanmean(offset_probs);
	end
	
	
	snps.Signature = snps.Signature(1:found);
	snps.Score = snps.Score(1:found);
	snps.SupportingReads = snps.SupportingReads(1:found);
	snps.TotalReads = snps.TotalReads(1:found);
	
	mutations.SNPs{s} = snps;
end

