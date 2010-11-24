function ensembl = read_ensembl_gtf(filepath)

if nargin == 0
	files = dir();
	for k = 1:length(files)
		if regexpi(files(k).name, '.*\.\d+\.gtf')
			filepath = files(k).name;
			break;
		end
	end
end

gene_prealloc = 100000;
tx_prealloc = 200000;

N_gene = 0;
N_tx = 0;

ensembl = struct;
ensembl.Genes = struct;
ensembl.Genes.ID = cell(gene_prealloc, 1);
ensembl.Genes.Name = cell(gene_prealloc, 1);
ensembl.Transcripts = struct;
ensembl.Transcripts.ID = cell(tx_prealloc, 1);
ensembl.Transcripts.Name = cell(tx_prealloc, 1);
ensembl.Transcripts.Exons = cell(tx_prealloc, 1);
ensembl.Transcripts.Gene = nan(tx_prealloc, 1);
ensembl.Transcripts.Chromosome = nan(tx_prealloc, 1);
ensembl.Transcripts.Strand = repmat(' ', tx_prealloc, 1);

gene_map = containers.Map;
transcript_map = containers.Map;

fid = fopen(filepath);
data = textscan(fid, '%s %s %s %d %d %*s %s %*s %s', 'Delimiter', '\t');
fclose(fid);

chromosome = data{1};
source = data{2};
type = data{3};
left = data{4};
right = data{5};
strand = data{6};
annot = data{7};

progress = Progress;

for k = 1:length(chromosome)
	chr = chromosome_sym2num(chromosome{k});
	if isnan(chr), continue, end
		
	if ~strcmp(type{k}, 'exon'), continue, end
	
	tokens = regexp(annot{k}, ['gene_id "(.+?)"; transcript_id "(.+?)";.*' ...
		'exon_number "(\d+?)"; gene_name "(.+?)"; transcript_name "(.+?)";'],...
		'tokens');
	if length(tokens) ~= 1
		annot{k}
		continue;
	end
	
		token = tokens{1};
		gene_id = token{1};
		tx_id = token{2};
		exon_num = str2double(token{3});
		gene_name = token{4};
		tx_name = token{5};
		
		if ~gene_map.isKey(gene_id)
			N_gene = N_gene + 1;
			gene_idx = N_gene;
			ensembl.Genes.ID{gene_idx} = gene_id;
			ensembl.Genes.Name{gene_idx} = gene_name;
			gene_map(gene_id) = gene_idx;
		else
			gene_idx = gene_map(gene_id);
		end
		
		if ~transcript_map.isKey(tx_id)
			N_tx = N_tx + 1;
			tx_idx = N_tx;
			ensembl.Transcripts.ID{tx_idx} = tx_id;
			ensembl.Transcripts.Name{tx_idx} = tx_name;
			ensembl.Transcripts.Gene(tx_idx) = gene_idx;
			ensembl.Transcripts.Exons{tx_idx} = [NaN NaN];
			ensembl.Transcripts.Chromosome(tx_idx) = chr;
			ensembl.Transcripts.Strand(tx_idx, 1) = strand{k};
			transcript_map(tx_id) = tx_idx;
		else
			tx_idx = transcript_map(tx_id);
		end
		
		ensembl.Transcripts.Exons{tx_idx}(exon_num, :) = [left(k) right(k)];
		
	progress.update(k / length(chromosome));
end

N_gene
N_tx

ensembl.Genes.ID = ensembl.Genes.ID(1:N_gene);
ensembl.Genes.Name = ensembl.Genes.Name(1:N_gene);
ensembl.Transcripts.ID = ensembl.Transcripts.ID(1:N_tx);
ensembl.Transcripts.Name = ensembl.Transcripts.Name(1:N_tx);
ensembl.Transcripts.Exons = ensembl.Transcripts.Exons(1:N_tx);
ensembl.Transcripts.Gene = ensembl.Transcripts.Gene(1:N_tx);
ensembl.Transcripts.Chromosome = ensembl.Transcripts.Chromosome(1:N_tx);
ensembl.Transcripts.Strand = ensembl.Transcripts.Strand(1:N_tx);

