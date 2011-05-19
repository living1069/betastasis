function gtf = read_ucsc_exon_gtf_new(filepath)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

fid = fopen(filepath);
data = textscan(fid, '%s %*s %s %d %d %*s %s %*s %s', 'Delimiter', '\t');
fclose(fid);

chromosome = chromosome_sym2num(data{1});
type = data{2};
position = [data{3} data{4}];
strand = data{5};
info = data{6};

progress = Progress;

gtf.Exons = cell(length(transcripts.Name), 1);
gtf.Chromosome = nan(length(transcripts.Name), 1);
gtf.Strand = repmat(' ', length(transcripts.Name), 1);

for f = 1:length(info)
	if ~strcmp(type{f}, 'exon'), continue, end
	
	tokens = regexpi(info{f}, 'transcript_id "(.+?)"', 'tokens');
	if length(tokens) ~= 1, fprintf(1, 'strange...\n'); continue, end
	
	token = tokens{1}; tx_idx = transcript_idx(token{1});
	if isnan(tx_idx)
		%fprintf(1, 'Unknown transcript %s.\n', token{1});
		continue;
	end
	
	gtf.Exons{tx_idx}(end+1, :) = position(f, :);
	if strand{f} == '-'
		gtf.Exons{tx_idx} = gtf.Exons{tx_idx}(end:-1:1, :);
	end

	gtf.Strand(tx_idx) = strand{f};
	gtf.Chromosome(tx_idx) = chromosome(f);
	
	progress.update(f / length(info));
end



% Use this code to augment the organism data structure with exon coordinates.
if 0 
	exons.Position = nan(length(exons.ID), 2);
	genes.Chromosome = nan(length(genes.Name), 1);
	genes.Strand = repmat(' ', length(genes.Name), 1);
	genes.Position = nan(length(genes.Name), 2);
	
	for tx = 1:length(organism.Transcripts.Exons)
		
		g = organism.Transcripts.Gene(tx);
		genes.Chromosome(g) = gtf.Chromosome(tx);
		genes.Strand(g) = gtf.Strand(tx);
		
		ex = organism.Transcripts.Exons{tx};
		if length(ex) == 0 || size(gtf.Exons{tx}, 1) ~= length(ex)
			continue;
		end
		
		exons.Position(ex, :) = gtf.Exons{tx};
		
		genes.Position(g, 1) = ...
			min(genes.Position(g, 1), min(gtf.Exons{tx}(:, 1)));
		genes.Position(g, 2) = ...
			max(genes.Position(g, 2), max(gtf.Exons{tx}(:, 2)));
	end
end

