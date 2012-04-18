function mutations = mirna_mutations(reads)

global organism;
mirnas = organism.miRNA;

threshold = 0.05;

mirna_len = zeros(length(mirnas.Sequence), 1);
for k = 1:length(mirna_len)
	mirna_len(k) = length(mirnas.Sequence{k});
end

S = length(reads.Raw);

mutations = struct;
mutations.Mutations = cell(1, S);

mutations.Meta = reads.Meta;
mutations.Meta.Type = 'miRNA mutations';
mutations.Meta.Organism = organism.Name;

isoform_range = min(mirna_len):max(mirna_len);

for s = 1:S
	extracted = extract_reads(filter_query(reads, s));
	nuc_reads = cell(length(mirnas.Name), 1);
	
	for iso_len = isoform_range
		
		% Build a Bowtie index for all isoforms of this length.
		index_mirnas = [];
		index_seqs = {};
		
		for m = 1:length(mirna_len)
			if iso_len ~= mirna_len(m), continue, end
			
			seq = mirnas.Sequence{m};
			
			index_mirnas(end+1, 1) = m;
			index_seqs{end+1, 1} = mirnas.Sequence{m};
		end
		
		fprintf(1, 'Trimming reads to a length of %d bases...\n', iso_len);
		trimmed = trim_reads(extracted, iso_len);
		
		color = ~isempty(regexpi(trimmed.Meta.Sequence.Space{1}, 'color'));
		
		al = align_reads(trimmed, index_seqs, 'MaxMismatches', 1+1*color, ...
			'Columns', 'target,offset,sequence', 'AllowAlignments', 1);
		
		targets = index_mirnas(str2double(al.Target));
		
		fprintf(1, 'Updating nucleotide frequencies...\n');
		for k = 1:length(targets)
			m = targets(k);
			
			if isempty(nuc_reads{m})
				nuc_reads{m} = zeros(4, length(mirnas.Sequence{m}), 'single');
			end

			% Linear indexing is used to efficiently increment the matrix
			% elements that correspond to the read sequence bases.
			pos = (al.Offset(k):al.Offset(k) + length(al.Sequence{k}) - 1) - 1;
			valid = al.Sequence{k} ~= 'N';   % Don't count unknown bases.
			lidx = pos(valid) * 4 + int32(nt2int(al.Sequence{k}(valid)));
			nuc_reads{m}(lidx) = nuc_reads{m}(lidx) + 1;
		end
	end
	
	fprintf(1, 'Looking for mutations...\n');
	
	% We preallocate the vectors to avoid triggering zillions of dynamic
	% memory allocations during the SNP discovery process.
	mut = struct;
	mut.Signature = cell(1000, 1);
	mut.Score = zeros(1000, 1);
	mut.SupportingReads = zeros(1000, 1);
	mut.TotalReads = zeros(1000, 1);
	
	progress = Progress;
	found = 0;
	
	for k = 1:length(mirnas.Name)
		progress.update(k / length(mirnas.Name));

		pwm = nuc_reads{k};
		if isempty(pwm), continue, end
		
		totals = sum(pwm, 1);
		seq = mirnas.Sequence{k};
		nuc_nums = nt2int(seq);
		
		scores = tanh(pwm / 100) .* pwm ./ repmat(totals, 4, 1);
		
		for p = 1:length(seq)
			for b = setdiff(1:4, nuc_nums(p))
				%if pwm(b, p) < 10, continue, end

				if scores(b, p) > threshold
					found = found + 1;
					mut.Signature{found} = sprintf('%s:%d:%s>%s', ...
						mirnas.Name{k}, p, seq(p), int2nt(b));
					mut.Score(found) = scores(b, p);
					mut.SupportingReads(found) = pwm(b, p);
					mut.TotalReads(found) = totals(p);
				end
			end
		end
	end
	
	mut.Signature = mut.Signature(1:found);
	mut.Score = mut.Score(1:found);
	mut.SupportingReads = mut.SupportingReads(1:found);
	mut.TotalReads = mut.TotalReads(1:found);
	
	mutations.Mutations{s} = mut;
end

