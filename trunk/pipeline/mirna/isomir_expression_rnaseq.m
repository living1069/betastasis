function isomirs = isomir_expression_rnaseq(reads, adapters)

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

S = length(reads.Raw);

isomirs = struct;
isomirs.Isoforms = cell(length(mirnas.Name), S);

isomirs.Meta = reads.Meta;
isomirs.Meta.Type = 'miRNA isoforms';
isomirs.Meta.Organism = organism.Name;

num_isoforms = 3 + 1 + 4 + 16 + 64;

for s = 1:S
	for m = 1:length(mirnas.Name)
		isomirs.Isoforms{m, s} = struct;
		isomirs.Isoforms{m, s}.Sequence = {};
		isomirs.Isoforms{m, s}.Count = zeros(num_isoforms, 1);
	end
end

isoform_range = min(mirna_len)-3:28;

for s = 1:S
	extracted = extract_reads(filter_query(reads, s));
	
	for iso_len = isoform_range
		
		% Build a Bowtie index for all isoforms of this length.
		isoforms = zeros(0, 2);
		isoform_seqs = {};
		
		for m = 1:length(mirna_len)
			if iso_len < mirna_len(m)-3 || iso_len > mirna_len(m)+3
				continue;
			end
			
			seq = mirnas.Sequence{m};
			
			if iso_len <= mirna_len(m)
				isomirs.Isoforms{m, s}.Sequence{end+1, 1} = seq(1:iso_len);
				isoform_seqs{end+1, 1} = seq(1:iso_len);
				isoforms(end+1, :) = ...
					[m, length(isomirs.Isoforms{m, s}.Sequence)];
			else
				kmers = dec2base(0:4^(iso_len - mirna_len(m))-1, 4);
				kmers(kmers == '0') = 'A';
				kmers(kmers == '1') = 'C';
				kmers(kmers == '2') = 'G';
				kmers(kmers == '3') = 'T';
				
				for k = 1:size(kmers, 1)
					isomirs.Isoforms{m, s}.Sequence{end+1,1} = [seq kmers(k,:)];
					isoform_seqs{end+1, 1} = [seq kmers(k,:)];
					isoforms(end+1, :) = ...
						[m, length(isomirs.Isoforms{m, s}.Sequence)];
				end
			end

			% Suffix the reference sequences with 3' adapters before alignment.
			isoforms_a = repmat(isoforms, length(adapters), 1);
			isoform_seqs_a = cell(length(adapters) * length(isoform_seqs), 1);
			
			for a = 1:length(adapters)
				isoform_seqs_a( ...
					(a-1)*length(isoform_seqs)+1:a*length(isoform_seqs)) = ...
					strcat(isoform_seqs, adapters{a});
			end

			%isoform_seq = seq(1:min(iso_len, mirna_len(m))); 
			%if iso_len > mirna_len(m)
			%	[i, j] = find(pre_mirnas.Matures == m);
			%	
			%	% If the mature microRNA can originate from many different
			%	% pre-miRNA, we don't try to look for longer isoforms.
			%	if length(i) ~= 1, continue; end
			%	
			%	preseq = pre_mirnas.Sequence{i};
			%	offset = pre_mirnas.MatureOffsets(i, j);
			%	
			%	mature_range = offset+mirna_len(m): ...
			%		offset+mirna_len(m)+(iso_len-mirna_len(m))-1;
			%	
			%	if max(mature_range) > length(preseq), continue; end
			%	
			%	isoform_seq = [isoform_seq preseq(mature_range)];
			%end
		end
		
		fprintf(1, 'Trimming reads to a length of %d + %d bases...\n', ...
			iso_len, adapter_len);
		trimmed = trim_reads(extracted, iso_len + adapter_len);
		
		al = align_reads(trimmed, isoform_seqs_a, 'MaxMismatches', 0, ...
			'Columns', 'target,offset,sequence', 'AllowAlignments', 1);
		
		iso_index = isoforms_a(str2double(al.Target), :);
		
		fprintf(1, 'Updating isoform counts...\n');
		for k = 1:size(iso_index, 1)
			m = iso_index(k, 1);
			l = iso_index(k, 2);
			isomirs.Isoforms{m, s}.Count(l) = ...
				isomirs.Isoforms{m, s}.Count(l) + 1;
		end
	end
end

