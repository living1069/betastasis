function mutations = mirna_mutations(reads, adapter)

global organism;
mirnas = organism.miRNA;

threshold = 0.05;

mirna_seq = mirnas.Sequence;
mirna_len = zeros(length(mirnas.Sequence), 1);
for k = 1:length(mirna_len)
	mirna_len(k) = length(mirna_seq{k});
end

adapter_len = length(adapter);
if adapter_len < 2
	error(['At least the first two nucleotides of the adapter sequence ' ...
	       'must be known.']);
end

shorter_isoforms = 3;
longer_isoforms = 3;

seq_files = seq_resource_files(reads);

% Figure out what kind of data we're dealing with.
[color, quality] = seq_read_type(seq_files{1});
if quality
	error 'FASTQ files not supported yet.';
end

max_mismatches = 1;
if color, max_mismatches = max_mismatches + 1; end

color_option = '';
index_suffix = '';
if color
	color_option = '-C';
	index_suffix = '_colorspace';
end


mutations = struct;
mutations.SNPs = cell(1, length(seq_files));
mutations.Isoforms = cell(length(mirnas.Name), length(seq_files));

for s = 1:length(seq_files)
	for m = 1:length(mirnas.Name)
		mutations.Isoforms{m, s} = struct;
		mutations.Isoforms{m, s}.Start = ...
			ones(shorter_isoforms + longer_isoforms + 1, 1);
		mutations.Isoforms{m, s}.End = ...
			mirna_len(m)-shorter_isoforms:mirna_len(m)+longer_isoforms;
		mutations.Isoforms{m, s}.Count = ...
			zeros(shorter_isoforms + longer_isoforms + 1, 1);
	end
end

mutations.Meta = struct;
if isstruct(reads), mutations.Meta = reads.Meta; end

mutations.Meta.Type = 'miRNA mutations';
mutations.Meta.Organism = organism.Name;
mutations.Meta.miRNAVersion = organism.miRNAVersion;

isoform_range = min(mirna_len)-shorter_isoforms:max(mirna_len)+longer_isoforms;

fasta_tmp = ptemp();
index_tmp = ptemp();
trimmed_reads_tmp = ptemp();

for seq_file = 1:length(seq_files)
	
	nuc_reads = cell(length(mirnas.Name), 1);
	
	for iso_len = isoform_range
		
		% Build a Bowtie index for all isoforms of this length.
		fasta_fid = fopen(fasta_tmp, 'w');
		
		isoform_map = containers.Map;
		
		for k = 1:length(mirna_len)
			if ~ismember(iso_len, ...
				mirna_len(k)-shorter_isoforms:mirna_len(k)+longer_isoforms)
				continue;
			end
			
			seq = mirna_seq{k};
			isoform_seq = seq(1:min(iso_len, mirna_len(k))); 
			if iso_len > mirna_len(k)
				[i, j] = find(organism.pre_miRNA.Matures == k);
				
				% If the mature microRNA can originate from many different
				% pre-miRNA, we don't try to look for longer isoforms.
				if length(i) ~= 1
					mutations.Isoforms{k, seq_file}.Count( ...
						iso_len - mirna_len(k) + shorter_isoforms + 1) = NaN;
					continue;
				end
				
				preseq = organism.pre_miRNA.Sequence{i};
				offset = organism.pre_miRNA.MatureOffsets(i, j);
				
				mature_range = offset+mirna_len(k): ...
					offset+mirna_len(k)+(iso_len-mirna_len(k))-1;
				
				if max(mature_range) > length(preseq)
					mutations.Isoforms{k, seq_file}.Count( ...
						iso_len - mirna_len(k) + shorter_isoforms + 1) = NaN;
					continue;
				end
					
				isoform_seq = [isoform_seq preseq(mature_range)];
			end
			isoform_seq = [isoform_seq adapter];
			
			name = sprintf('%s:%d', mirnas.Name{k}, iso_len);
			fprintf(fasta_fid, '>%s\n%s\n', name, isoform_seq);
			
			isoform_map(name) = ...
				[k, iso_len - mirna_len(k) + shorter_isoforms + 1];
		end
		
		fclose(fasta_fid);

		[status, ~] = unix(sprintf('%s/tools/bowtie/bowtie-build %s %s %s%s',...
			ppath, color_option, fasta_tmp, index_tmp, index_suffix));

		if status ~= 0, error 'Bowtie index construction failed.'; end
		system(['rm ' fasta_tmp]);
		
		fprintf(1, 'Trimming reads to a length of %d + %d bp...\n', ...
			iso_len, adapter_len);
		[status, ~] = unix(sprintf( ...
			'%s/sources/sequencing/trim_reads.py %s %d > %s', ...
			ppath, seq_files{seq_file}, ...
			iso_len + adapter_len + 1 * color, trimmed_reads_tmp));

		fprintf(1, 'Aligning reads using Bowtie...\n');
		[alignments_tmp, ~] = bowtie_align(trimmed_reads_tmp, index_tmp, ...
			sprintf('-p4 -v%d -B1 --suppress 1,2,6,7,8', max_mismatches));
			
		system(['rm ' trimmed_reads_tmp]);
		system(['rm ' index_tmp '*']);

		alignments_file = fopen(alignments_tmp);
		data = textscan(alignments_file, '%s %d %s');
		fclose(alignments_file);
		
		isoform_names = data{1};
		read_offsets = data{2};
		read_sequences = data{3};
		
		isoform_indices = isoform_map.values(isoform_names);
		
		fprintf(1, 'Updating nucleotide frequencies...\n');
		for k = 1:length(isoform_indices)
			idx = isoform_indices{k}; idx = idx(1);
			
			if isempty(nuc_reads{idx})
				nuc_reads{idx} = zeros(4, length(mirna_seq{idx}) + ...
					adapter_len + 20, 'single');
			end

			% Linear indexing is used to efficiently increment the matrix
			% elements that correspond to the read sequence bases.
			pos = ((read_offsets(k): ...
				read_offsets(k) + length(read_sequences{k}) - 1) - 1) * 4 + ...
				int32(nt2int(read_sequences{k}));
			nuc_reads{idx}(pos) = nuc_reads{idx}(pos) + 1;
		end

		fprintf(1, 'Updating isoform counts...\n');
		for k = 1:length(isoform_indices)
			idx = isoform_indices{k};
			mutations.Isoforms{idx(1), seq_file}.Count(idx(2)) = ...
				mutations.Isoforms{idx(1), seq_file}.Count(idx(2)) + 1;
		end
	end
	
	fprintf(1, 'Looking for SNPs...\n');
	
	% We preallocate the vectors to avoid triggering zillions of dynamic
	% memory allocations during the SNP discovery process.
	snps = struct('Signature', { cell(1000, 1) }, ...
	              'Score', zeros(1000, 1), ...
				  'SupportingReads', zeros(1000, 1), ...
				  'TotalReads', zeros(1000, 1));
	
	progress = Progress;
	found = 0;
	
	for k = 1:length(mirnas.Name)
		progress.update(k / length(mirnas.Name));

		reads = nuc_reads{k};
		if isempty(reads), continue, end
		
		totals = sum(reads, 1);
		seq = mirna_seq{k};
		nuc_nums = nt2int(seq);
		
		scores = tanh(reads / 100) .* reads ./ repmat(totals, 4, 1);
		
		for s = 1:length(seq)
			for b = setdiff(1:4, nuc_nums(s))
				%if reads(b, s) < 10, continue, end

				if scores(b, s) > threshold
					found = found + 1;
					snps.Signature{found} = sprintf('%s:%d:%s>%s', ...
						mirnas.Name{k}, s, seq(s), int2nt(b));
					snps.Score(found) = scores(b, s);
					snps.SupportingReads(found) = reads(b, s);
					snps.TotalReads(found) = totals(s);
				end
			end
		end
	end
	
	snps.Signature = snps.Signature(1:found);
	snps.Score = snps.Score(1:found);
	snps.SupportingReads = snps.SupportingReads(1:found);
	snps.TotalReads = snps.TotalReads(1:found);
	
	mutations.SNPs{seq_file} = snps;
end
	
return;



