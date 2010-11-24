function nucleotide = colorspace_to_nucleotide(first_nuc, color_seq)

% Rows: ACGT. Columns: BGYR.
decode_matrix = [ 'ACGT';
                  'CATG';
				  'GTAC';
				  'TGCA' ];

nucleotide = first_nuc;

if color_seq(1) == 'A' || color_seq(1) == 'C' || color_seq(1) == 'G' || ...
   color_seq(1) == 'T'
    color_seq = color_seq(2:end);
end

color_seq = color_seq(3:end);
				
for k = 1:length(color_seq)
	if color_seq(k) == '.', break, end
	color_num = str2num(color_seq(k)) + 1;
	nucleotide(end + 1) = decode_matrix(nuc_to_idx(nucleotide(end)), color_num);
end

return;




function idx = nuc_to_idx(nuc)
nucleotides = 'ACGT';
idx = strfind(nucleotides, nuc);
return;

