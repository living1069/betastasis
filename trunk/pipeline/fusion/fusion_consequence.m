function [] = fusion_consequence(exon_5p, exon_3p)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

tokens = regexpi(exon_5p, '(.+)\[(.+)\]', 'tokens');
if length(tokens) ~= 1, error 'Invalid 5'' exon name.'; end
token = tokens{1}; gene_5p = token{1}; exon_5p = token{2};
	
tokens = regexpi(exon_3p, '(.+)\[(.+)\]', 'tokens');
if length(tokens) ~= 1, error 'Invalid 3'' exon name.'; end
token = tokens{1}; gene_3p = token{1}; exon_3p = token{2};

left_gene = gene_idx(gene_5p);
right_gene = gene_idx(gene_3p);

if isnan(left_gene), error 'Unknown 5'' gene.'; end
if isnan(right_gene), error 'Unknown 3'' gene.'; end

exons_5p = [];
for tx = genes.Transcripts(left_gene, 1:genes.TranscriptCount(left_gene))
	exons_5p = [exons_5p; transcripts.Exons{tx}];
end
exons_3p = [];
for tx = genes.Transcripts(right_gene, 1:genes.TranscriptCount(right_gene))
	exons_3p = [exons_3p; transcripts.Exons{tx}];
end

exons_5p = unique(exons_5p);
exons_3p = unique(exons_3p);

left_exon = exons_5p(strcmp(exon_5p, exons.ID(exons_5p)));
right_exon = exons_3p(strcmp(exon_3p, exons.ID(exons_3p)));



roles_known = false;
left_exon_frames = [];
right_exon_frames = [];

left_exon_roles = {};
left_txs = genes.Transcripts(left_gene, ...
	1:genes.TranscriptCount(left_gene));
for t = 1:length(left_txs)
	cds = transcripts.CDS(left_txs(t), :);
	if all(isnan(cds))
		left_exon_roles{end+1, 1} = 'Non-coding';
		continue;
	end
	
	ex = find(transcripts.Exons{left_txs(t)} == left_exon);
	if isempty(ex), continue, end
		
	exon_pos = transcripts.ExonPos{left_txs(t)};
	exon_pos = exon_pos(ex, :);
	
	if exon_pos(2) < cds(1)
		left_exon_roles{end+1, 1} = '5'' UTR';
	elseif exon_pos(2) >= cds(1) && exon_pos(2) <= cds(2)
		left_exon_roles{end+1, 1} = 'CDS';
	elseif exon_pos(2) > cds(2)
		left_exon_roles{end+1, 1} = '3'' UTR';
	end
	
	left_exon_frames(end+1) = mod(exon_pos(2) - cds(1) + 1, 3);
end


right_exon_roles = {};
right_txs = genes.Transcripts(right_gene, ...
	1:genes.TranscriptCount(right_gene));
for t = 1:length(right_txs)
	cds = transcripts.CDS(right_txs(t), :);
	if all(isnan(cds))
		right_exon_roles{end+1, 1} = 'Non-coding';
		continue;
	end
	
	ex = find(transcripts.Exons{right_txs(t)} == right_exon);
	if isempty(ex), continue, end
		
	exon_pos = transcripts.ExonPos{right_txs(t)};
	exon_pos = exon_pos(ex, :);
	
	if exon_pos(2) < cds(1)
		right_exon_roles{end+1, 1} = '5'' UTR';
	elseif exon_pos(2) >= cds(1) && exon_pos(2) <= cds(2)
		right_exon_roles{end+1, 1} = 'CDS';
	elseif exon_pos(2) > cds(2)
		right_exon_roles{end+1, 1} = '3'' UTR';
	end
	
	right_exon_frames(end+1) = mod(exon_pos(1) - cds(1), 3);
end


% Now we check that we have a consensus on the roles that the left and
% right exon play in the RNA transcripts.
left_exon_roles = unique(left_exon_roles);
right_exon_roles = unique(right_exon_roles);
if length(left_exon_roles) == 1 && length(right_exon_roles) == 1
	roles_known = true;
	left_exon_role = left_exon_roles{1};
	right_exon_role = right_exon_roles{1};
else
	fprintf(1, 'No agreement on exon roles for %s:%s.\n', ...
		genes.Name{left_gene}, genes.Name{right_gene});
	left_exon_roles
	right_exon_roles
end

left_exon_frames = unique(left_exon_frames);
right_exon_frames = unique(right_exon_frames);

if roles_known
	if strcmpi(left_exon_role, 'non-coding') && ...
		strcmpi(right_exon_role, 'non-coding')
		fprintf(1, 'Both genes are non-coding.\n');
	
	elseif strcmpi(left_exon_role, 'non-coding')
		fprintf(1, 'Non-coding gene %s is fused to the %s of %s.\n', ...
			genes.Name{left_gene}, right_exon_role, genes.Name{right_gene});
	elseif strcmpi(right_exon_role, 'non-coding')
		fprintf(1, '%s of %s is fused to non-coding gene %s.\n', ...
			left_exon_role, genes.Name{left_gene}, genes.Name{right_gene});
	else
		fprintf(1, '%s of %s is fused to the %s of %s.\n', ...
			left_exon_role, genes.Name{left_gene}, ...
			right_exon_role, genes.Name{right_gene});
	end
	
	if strcmpi(left_exon_role, 'cds') && strcmpi(right_exon_role, 'cds')
		if length(left_exon_frames) > 1 || length(right_exon_frames) > 1
			fprintf(1, 'No consensus on frameshift for %s:%s.\n', ...
				genes.Name{left_gene}, genes.Name{right_gene});
		elseif left_exon_frames == right_exon_frames
			fprintf(1, 'No frameshift.\n');
		else
			fprintf(1, ['Fusion causes a frameshift.\nLeft exon terminates ' ...
				'at codon boundary +%d.\nRight exon begins at codon ' ...
				'boundary +%d.\n'], left_exon_frames, right_exon_frames);
				
			% OK, we now know that a frameshift is caused. Next we should check
			% if a new stop codon is introduced somewhere downstream.
			for tx = right_txs'
				if any(isnan(transcripts.CDS(tx, :))), continue, end
				
				ex_idx = find(transcripts.Exons{tx} == right_exon);
				if isempty(ex_idx), continue, end
					
				ex_start = transcripts.ExonPos{tx}(ex_idx, 1);
				left_exon_seq = exons.Sequence{left_exon};
				
				shift_seq = [left_exon_seq(end-left_exon_frames+1:end), ...
					transcripts.Sequence{tx}(ex_start:end)];
				
				shift_aa = nt2aa(shift_seq, 'AlternativeStartCodons', false);
				stops = find(shift_aa == '*');
				
				fprintf(1, 'Normally stop codon after %d aa.\n', ...
					ceil((transcripts.CDS(tx, 2) - ex_start + 1) / 3));
				
				if isempty(stops)
					fprintf(1, 'No stop codon if 3'' side is from %s.\n', ...
						transcripts.Name{tx});
				else
					fprintf(1, ...
						'Stop codon after %d aa if 3'' side from %s:\n%s\n', ...
						stops(1), transcripts.Name{tx}, shift_aa);
				end
			end
			
		end
	end
end


