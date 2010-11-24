function [] = nucleotide_to_colorspace(seq)

% Rows: ACGT. Columns: BGYR.
decode_matrix = [ 'ACGT';
                  'CATG';
				  'GTAC';
				  'TGCA' ];
				  
seq = upper(seq);
  
% This could use some optimization if it ever sees more use.
for pre_nuc = 1:4
    prev_nuc = pre_nuc;
    for k = 1:length(seq)
        cseq(k) = find(decode_matrix(prev_nuc, :) == seq(k)) - 1;
        prev_nuc = nuc_to_idx(seq(k));
    end
    
    fprintf(1, '%s%s\n', decode_matrix(1, pre_nuc), char(cseq + '0'));
    %fprintf(1, '%s\n', char(cseq + '0'));
end


function idx = nuc_to_idx(nuc)
nucleotides = 'ACGT';
idx = strfind(nucleotides, nuc);
return;


