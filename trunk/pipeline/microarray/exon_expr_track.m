
% Author: Matti Annala <matti.annala@tut.fi>

function [] = exon_expr_track(exon_expr, test_gene, ref_exon, ref_gene, ...
	track_dir, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;


test_exon_log_expr = log2(test_exon.Mean);
test_gene_log_expr = log2(test_gene.Mean);
ref_exon_log_expr = log2(ref_exon.Mean);
ref_gene_log_expr = log2(ref_gene.Mean);

if size(test_exon_log_expr, 2) ~= size(test_gene_log_expr, 2) || ...
	size(ref_exon_log_expr, 2) ~= size(ref_gene_log_expr, 2)
	error 'Exon and gene expression matrices must have equal dimensions.';
end

test_splice_index = nan(size(test_exon_log_expr));
ref_splice_index = nan(size(ref_exon_log_expr));

E = size(test_exon_log_expr, 1);

for k = 1:E
	g = exons.Gene(k);
	test_splice_index(k, :) = ...
		test_exon_log_expr(k, :) - test_gene_log_expr(g, :);
	ref_splice_index(k, :) = ...
		ref_exon_log_expr(k, :) - ref_gene_log_expr(g, :);
end



fprintf(1, 'Writing the exon expression tracks to %s...\n', track_dir);
[~, ~] = mkdir(track_dir);
[~, ~] = mkdir([track_dir '/splice_index']);

fid = fopen([track_dir '/genelist.json'], 'w');
fprintf(fid, '{ "genes": [ ');

for g = 1:length(genes.Name)
	if g == 1
		fprintf(fid, '"%s"', genes.Name{g});
	else
		fprintf(fid, ', "%s"', genes.Name{g});
	end
end

fprintf(fid, ' ] }\n');
fclose(fid);



for g = 1:length(genes.Name)
	[~, ~] = mkdir([track_dir '/splice_index/' lower(genes.Name{g}(1))]);
	fid = fopen([track_dir '/splice_index/' lower(genes.Name{g}(1)) '/' ...
		genes.Name{g} '.json'], 'W');
	
	fprintf(fid, '{\n');
		
	tx_indices = find(transcripts.Gene == g);
	gene_exons = [];

	for t = 1:length(tx_indices)
		tx_idx = tx_indices(t);
		tx_seq = transcripts.Sequence{tx_idx};
		
		fprintf(fid, '"%s": {\n', transcripts.Name{tx_idx});
		fprintf(fid, '\t"length": %d,\n', length(tx_seq));
		fprintf(fid, '\t"exons": [');
		
		tx_exons = transcripts.Exons{tx_idx};
		for m = 1:length(tx_exons)
			exon_id = exons.ID{tx_exons(m)};
			if m == length(tx_exons)
				fprintf(fid, '"%s"],\n', exon_id);
			else
				fprintf(fid, '"%s", ', exon_id);
			end
		end
		
		if length(tx_exons) == 0
			fprintf(fid, '],\n');
		end
		
		fprintf(fid, '\t"exon_pos": [');
		
		tx_exon_pos = transcripts.ExonPos{tx_idx};
		for m = 1:size(tx_exon_pos, 1)
			left = tx_exon_pos(m, 1);
			right = tx_exon_pos(m, 2);
			
			if m == size(tx_exon_pos, 1)
				fprintf(fid, '%d],\n', left);
			else
				fprintf(fid, '%d, ', left);
			end
		end
		
		if size(tx_exon_pos, 1) == 0
			fprintf(fid, '1],\n');
		end
		
		cds = transcripts.CDS(tx_idx, :);
		if ~any(isnan(cds))
			cds(isnan(cds)) = 0;
			fprintf(fid, '\t"cds": [%d, %d],\n', cds(1), cds(2));
		end
		
		fprintf(fid, '},\n');
		
		gene_exons = [gene_exons; tx_exons];
	end

	gene_exons = unique(gene_exons);
	
	fprintf(fid, '"exons": [');
	for m = 1:length(gene_exons)
		exon_id = exons.ID{gene_exons(m)};
		if m == length(gene_exons)
			fprintf(fid, '"%s"],\n', exon_id);
		else
			fprintf(fid, '"%s", ', exon_id);
		end
	end
	
	if length(gene_exons) == 0
		fprintf(fid, '],\n');
	end
	
	
	fprintf(fid, '"splice_index": [\n');
	for m = 1:length(gene_exons)
		ex = gene_exons(m);
		fprintf(fid, '\t[');
		fprintf(fid, '%f', 
		if m == length(gene_exons)
			fprintf(fid, '"%s"],\n', exon_id);
		else
			fprintf(fid, '"%s", ', exon_id);
		end
	end
	
	if length(gene_exons) == 0
		fprintf(fid, '],\n');
	end


		
	fprintf(fid, '}\n');
	fclose(fid);
end










