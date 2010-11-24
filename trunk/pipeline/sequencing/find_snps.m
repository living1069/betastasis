function mutations = find_snps(reads, transcripts, adapter)

global organism;

threshold = 0.05;

seq_files = seq_resource_files(reads);

if nargin < 2
	transcripts = organism.Transcripts;
end

if nargin < 3, adapter = ''; end
adapter = upper(adapter);
adapter_len = length(adapter);

fasta_tmp = ptemp();
index_tmp = ptemp();

output = fopen(fasta_tmp, 'W');

for k = 1:length(transcripts.Name)
	seq = transcripts.Sequence{k};
	fprintf(output, '>%s\n%s%s\n', transcripts.Name{k}, seq, adapter);
end

fclose(output);

[color, quality] = seq_read_type(seq_files{1});

color_option = '';
index_suffix = '';
if color
	color_option = '-C';
	index_suffix = '_colorspace';
end

[status, ~] = unix(sprintf('%s/tools/bowtie/bowtie-build %s %s %s%s', ...
	ppath, color_option, fasta_tmp, index_tmp, index_suffix));
	
system(['rm ' fasta_tmp]);

if status ~= 0, error 'Bowtie index construction failed.'; end

transcript_name_to_idx = containers.Map(transcripts.Name, ...
	num2cell(1:length(transcripts.Name)));
	
mutations = struct;
mutations.SNPs = {};

mutations.Meta = struct;
if isfield(reads, 'Meta'), mutations.Meta = reads.Meta; end
	
mutations.Meta.Type = 'SNPs';
mutations.Meta.Organism = organism.Name;
mutations.Meta.ScoreThreshold = repmat([threshold], length(seq_files), 1);

for seq_file = 1:length(seq_files)
	
	num_mismatches = 2;
	if color, num_mismatches = 2; end
	
	fprintf(1, 'Aligning reads to sequences + adapters...\n');
	[alignments_tmp, out] = bowtie_align(seq_files{seq_file}, ...
		index_tmp, sprintf('-p4 -v%d -m10 -B1 --norc --suppress 1,2,6,7,8', ...
		num_mismatches));

	fprintf(1, 'Reading aligned reads into memory...\n');
	alignments_file = fopen(alignments_tmp);
	data = textscan(alignments_file, '%s %d %s', 'Delimiter', '\t');
	fclose(alignments_file);

	%fprintf(1, 'Alignments stored at %s.\n', alignments_tmp);
	system(['rm ' alignments_tmp]);

	transcript_names = data{1};
	read_offsets = data{2};
	read_sequences = data{3};
	clear data;
	
	nuc_reads = cell(length(transcripts.Name), 1);
	
	fprintf(1, 'Calculating nucleotide frequencies...\n');
	transcript_indices = cell2mat( ...
		transcript_name_to_idx.values(transcript_names));
	
	for k = 1:length(transcript_names)
		m = transcript_indices(k);
		
		if isempty(nuc_reads{m})
			nuc_reads{m} = zeros(4, length(transcripts.Sequence{m}) + ...
				adapter_len + 20, 'single');
		end

		% Linear indexing is used to efficiently increment the matrix elements 
		% that correspond to the read sequence bases.
		pos = ((read_offsets(k): ...
			read_offsets(k) + length(read_sequences{k}) - 1) - 1) * 4 + ...
			int32(nt2int(read_sequences{k}));
		nuc_reads{m}(pos) = nuc_reads{m}(pos) + 1;
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
	end
	
	%offset_probs = offset_errors ./ offset_reads;
	%offset_probs = offset_probs / nanmean(offset_probs);
	
	fprintf(1, 'Looking for SNPs...\n');
	
	% We preallocate the vectors to avoid triggering zillions of dynamic
	% memory allocations during the SNP discovery process.
	snps = struct('Signature', { cell(10000, 1) }, ...
	              'Score', zeros(10000, 1), ...
				  'SupportingReads', zeros(10000, 1), ...
				  'TotalReads', zeros(10000, 1));
	
	progress = Progress;
	found = 0;
	
	for k = 1:length(transcripts.Name)
		progress.update(k / length(transcripts.Name));

		reads = nuc_reads{k};
		if isempty(reads), continue, end
		
		totals = sum(reads, 1);
		seq = transcripts.Sequence{k};
		nuc_nums = nt2int(seq);
		
		scores = tanh(reads / 100) .* reads ./ repmat(totals, 4, 1);
		
		for s = 1:length(seq)
			for b = setdiff(1:4, nuc_nums(s))
				%if reads(b, s) < 10, continue, end

				if scores(b, s) > threshold
					found = found + 1;
					snps.Signature{found} = sprintf('%s:%d:%s>%s', ...
						transcripts.Name{k}, s, seq(s), int2nt(b));
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

system(['rm -r ' index_tmp '*']);

