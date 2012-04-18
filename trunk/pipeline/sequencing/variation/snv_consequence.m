
function [consequences, flank_seq] = snv_consequence(snv, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

% Parse the variant 
tokens = regexpi(snv, '^chr(.+?):\s*(\d+)\s*:\s*(.+)>(.+)$', 'tokens');
if length(tokens) ~= 1, error 'Invalid SNV specification.'; end
	
token = tokens{1};
chr = chromosome_sym2num(token{1});
pos = str2double(token{2});
ref_allele = token{3};
alt_allele = token{4};

ref_allele_len = length(ref_allele);
alt_allele_len = length(alt_allele);




% Read 30 bases of flanking sequence from both sides of the variant site.
% The flanks are 30 + RAL bases long, where RAL = reference allele length.
% Both flanks incorporate the reference sequence.
left_flank = chromosomes.Sequence{chr}(pos-30:pos+ref_allele_len-1);
left_flank_rc = seqrcomplement(left_flank);
right_flank = chromosomes.Sequence{chr}(pos:pos+30+ref_allele_len-1);
right_flank_rc = seqrcomplement(right_flank);

flank_seq = lower([left_flank(1:30) right_flank]);
flank_seq(31:30+ref_allele_len) = upper(flank_seq(31:30+ref_allele_len));

% Sanity check that the reference allele actually matches with the genome build.
if right_flank(1:ref_allele_len) ~= ref_allele
	error 'Reference allele does not match with genome build.';
end





consequences = {};

overlap_genes = find(genes.Chromosome == chr & ...
	genes.Position(:, 1) - 10e3 <= pos & genes.Position(:, 2) + 10e3 >= pos)';

for g = overlap_genes
	txs = genes.Transcripts(g, 1:genes.TranscriptCount(g));
	gene_name = genes.Name{g};
	
	for tx = txs
		% We try to identify the site within transcript sequences, rather than
		% using annotated exon positions, because the latter may not be 100%
		% correct.
		
		tx_seq = upper(transcripts.Sequence{tx});
		
		if genes.Strand(g) == '+'
			offset = strfind(tx_seq, left_flank) + 30;
			if length(offset) ~= 1
				offset = strfind(tx_seq, right_flank);
			end
			tx_alleles = {ref_allele, alt_allele};
		elseif genes.Strand(g) == '-'
			offset = strfind(tx_seq, left_flank_rc);
			if length(offset) ~= 1
				offset = strfind(tx_seq, right_flank_rc) + 30;
			end
			tx_alleles = {seqcomplement(ref_allele), seqcomplement(alt_allele)};
		end
		
		if length(offset) ~= 1, continue; end
		
		tx_name = transcripts.Name{tx};
		
		if any(isnan(transcripts.CDS(tx, :)))
			consequences{end+1, 1} = ...
				sprintf('Exon of non-coding transcript %s', tx_name);
			continue;
		end
		
		cds_offset = offset - transcripts.CDS(tx, 1) + 1;
		cds_len = transcripts.CDS(tx, 2) - transcripts.CDS(tx, 1) + 1;
		
		if cds_offset < 1
			consequences{end+1, 1} = ...
				sprintf('5'' UTR of transcript %s (%s)', tx_name, gene_name);
		elseif cds_offset > cds_len
			consequences{end+1, 1} = ...
				sprintf('3'' UTR of transcript %s (%s)', tx_name, gene_name);
		elseif ref_allele_len == 1 && alt_allele_len == 1
			codon = floor((cds_offset - 1) / 3) + 1;
			ref_codon = tx_seq(transcripts.CDS(tx, 1) +((codon-1)*3:codon*3-1));
			ref_aa = nt2aa(ref_codon, 'AlternativeStartCodons', false);
			
			if ref_codon(mod(cds_offset - 1, 3) + 1) ~= tx_alleles{1}
				error 'Ref codon base does not match allele base.';
			end
			
			alt_codon = ref_codon;
			alt_codon(mod(cds_offset - 1, 3) + 1) = tx_alleles{2};
			alt_aa = nt2aa(alt_codon, 'AlternativeStartCodons', false);
			
			if ref_aa == alt_aa
				consequences{end+1, 1} = sprintf( ...
					['Synonymous mutation (%s -> %s) at codon %d ' ...
					'of transcript %s (%s)'], ref_codon, alt_codon, ...
					codon, tx_name, gene_name);
			elseif alt_aa == '*'
				consequences{end+1, 1} = sprintf( ...
					['Nonsense mutation (%s -> %s) at codon %d ' ...
					'of transcript %s (%s)'], ref_codon, alt_codon, ...
					codon, tx_name, gene_name);
			else
				consequences{end+1, 1} = sprintf( ...
					'Missense mutation %s%d%s in transcript %s (%s)', ...
					ref_aa, codon, alt_aa, tx_name, gene_name);
			end
		elseif mod(abs(alt_allele_len - ref_allele_len), 3) ~= 0
			% Indel within coding region of a gene, causes frameshift.
			consequences{end+1, 1} = sprintf( ...
				'Frameshift mutation in transcript %s (%s)', ...
				tx_name, gene_name);
		elseif alt_allele_len > ref_allele_len
			% Extra amino acids added.
			consequences{end+1, 1} = sprintf( ...
				'Insertion within CDS of transcript %s (%s)', ...
				tx_name, gene_name);
		elseif alt_allele_len < ref_allele_len
			% Amino acids deleted.
			consequences{end+1, 1} = sprintf( ...
				'Deletion within CDS of transcript %s (%s)', ...
				tx_name, gene_name);
		end
	end
end
	
