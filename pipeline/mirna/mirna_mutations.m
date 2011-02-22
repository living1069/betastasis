function mutations = mirna_mutations(reads, adapters)

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

threshold = 0.05;

mirna_len = zeros(length(mirnas.Sequence), 1);
for k = 1:length(mirna_len)
	mirna_len(k) = length(mirnas.Sequence{k});
end

if ischar(adapters), adapters = { adapters }; end
	
adapter_len = zeros(1, length(adapters));
for k = 1:length(adapter_len)
	adapter_len(k) = length(adapters{k});
end

if length(unique(adapter_len)) ~= 1
	error 'All adapter sequences must be of identical length.';
end
adapter_len = adapter_len(1);

shorter_isoforms = 3;
longer_isoforms = 3;

S = length(reads.Raw);

mutations = struct;
mutations.SNPs = cell(1, S);
mutations.Isoforms = cell(length(mirnas.Name), S);

mutations.Meta = reads.Meta;
mutations.Meta.Type = 'miRNA mutations';
mutations.Meta.Organism = organism.Name;

for s = 1:S
	for m = 1:length(mirnas.Name)
		mutations.Isoforms{m, s} = struct;
		mutations.Isoforms{m, s}.Start = ...
			ones(shorter_isoforms + longer_isoforms + 1, 1);
		mutations.Isoforms{m, s}.End = ...
			mirna_len(m)-shorter_isoforms:mirna_len(m)+longer_isoforms;
		mutations.Isoforms{m, s}.Count = ...
			nan(shorter_isoforms + longer_isoforms + 1, 1);
	end
end

isoform_range = min(mirna_len)-shorter_isoforms:max(mirna_len)+longer_isoforms;

for s = 1:S
	extracted = extract_reads(filter_query(reads, s));
	nuc_reads = cell(length(mirnas.Name), 1);
	
	for iso_len = isoform_range
		
		% Build a Bowtie index for all isoforms of this length.
		isoforms = zeros(0, 2);
		isoform_seqs = {};
		
		for m = 1:length(mirna_len)
			if iso_len < mirna_len(m)-shorter_isoforms || ...
				iso_len > mirna_len(m)+longer_isoforms
				continue;
			end
			
			seq = mirnas.Sequence{m};
			isoform_seq = seq(1:min(iso_len, mirna_len(m))); 
			if iso_len > mirna_len(m)
				[i, j] = find(pre_mirnas.Matures == m);
				
				% If the mature microRNA can originate from many different
				% pre-miRNA, we don't try to look for longer isoforms.
				if length(i) ~= 1, continue; end
				
				preseq = pre_mirnas.Sequence{i};
				offset = pre_mirnas.MatureOffsets(i, j);
				
				mature_range = offset+mirna_len(m): ...
					offset+mirna_len(m)+(iso_len-mirna_len(m))-1;
				
				if max(mature_range) > length(preseq), continue; end
					
				isoform_seq = [isoform_seq preseq(mature_range)];
			end
			
			for a = 1:length(adapters)
				isoform_seqs{end+1, 1} = [isoform_seq adapters{a}];
				isoforms(end+1, :) = ...
					[m, iso_len - mirna_len(m) + shorter_isoforms + 1];
			end
		end
		
		fprintf(1, 'Trimming reads to a length of %d + %d bases...\n', ...
			iso_len, adapter_len);
		trimmed = trim_reads(extracted, iso_len + adapter_len);
		
		al = align_reads(trimmed, isoform_seqs, 'MaxMismatches', 0, ...
			'Columns', 'target,offset,sequence', 'AllowAlignments', 1);
		
		iso_index = isoforms(str2double(al.Target), :);
		
		fprintf(1, 'Updating nucleotide frequencies...\n');
		for k = 1:size(iso_index, 1)
			m = iso_index(k, 1);
			if isempty(nuc_reads{m})
				nuc_reads{m} = zeros(4, length(mirnas.Sequence{m}) + ...
					adapter_len + 10, 'single');
			end

			% Linear indexing is used to efficiently increment the matrix
			% elements that correspond to the read sequence bases.
			pos = ((al.Offset(k): ...
				al.Offset(k) + length(al.Sequence{k}) - 1) - 1) * 4 + ...
				int32(nt2int(al.Sequence{k}));
			nuc_reads{m}(pos) = nuc_reads{m}(pos) + 1;
		end

		fprintf(1, 'Updating isoform counts...\n');
		for k = 1:size(iso_index, 1)
			m = iso_index(k, 1);
			l = iso_index(k, 2);
			mutations.Isoforms{m, s}.Count(l) = ...
				mutations.Isoforms{m, s}.Count(l) + 1;
		end
	end
	
	fprintf(1, 'Looking for SNPs...\n');
	
	% We preallocate the vectors to avoid triggering zillions of dynamic
	% memory allocations during the SNP discovery process.
	snps = struct;
	snps.Signature = cell(1000, 1);
	snps.Score = zeros(1000, 1);
	snps.SupportingReads = zeros(1000, 1);
	snps.TotalReads = zeros(1000, 1);
	
	progress = Progress;
	found = 0;
	
	for k = 1:length(mirnas.Name)
		progress.update(k / length(mirnas.Name));

		reads = nuc_reads{k};
		if isempty(reads), continue, end
		
		totals = sum(reads, 1);
		seq = mirnas.Sequence{k};
		nuc_nums = nt2int(seq);
		
		scores = tanh(reads / 100) .* reads ./ repmat(totals, 4, 1);
		
		for p = 1:length(seq)
			for b = setdiff(1:4, nuc_nums(p))
				%if reads(b, p) < 10, continue, end

				if scores(b, p) > threshold
					found = found + 1;
					snps.Signature{found} = sprintf('%s:%d:%s>%s', ...
						mirnas.Name{k}, p, seq(p), int2nt(b));
					snps.Score(found) = scores(b, p);
					snps.SupportingReads(found) = reads(b, p);
					snps.TotalReads(found) = totals(p);
				end
			end
		end
	end
	
	snps.Signature = snps.Signature(1:found);
	snps.Score = snps.Score(1:found);
	snps.SupportingReads = snps.SupportingReads(1:found);
	snps.TotalReads = snps.TotalReads(1:found);
	
	mutations.SNPs{s} = snps;
end

