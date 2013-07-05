


% BLADDER CANCER - TCGA BLCA
cd ~/datasets/tcga_blca/RNASeq/UNC__IlluminaHiSeq_RNASeq/Level_3
gene_expr = import_tcga_level3_rnaseq();
tumor = find(rx(gene_expr.meta.sample_id, '-01A-'));
normal = find(rx(gene_expr.meta.sample_id, '-(10|11)A-'));
gene_expr = filter(gene_expr, [tumor, normal]);
[~, order] = sort_nat(gene_expr.rows.gene);
gene_expr = filter_rows(gene_expr, order);
fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '~/rnaseq_gene_expression');

















% PROSTATE CANCER - PRENSNER (2011), CELL LINES
load ~/datasets/prensner_pca_rnaseq/gene_expression
%gene_expr.rows.gene = ensembl_to_hugo(gene_expr.rows.name, '~/organisms/homo_sapiens/ensembl_67/Homo_sapiens.GRCh37.67.gtf');
gene_expr.mean = quantilenorm(gene_expr.mean + 1);

valid = ~rx(gene_expr.rows.gene, '^(RP11|AC\d)');
gene_expr = filter_rows(gene_expr, valid);

count = tabulate(gene_expr.rows.gene);
valid = ismember(gene_expr.rows.gene, count(cell2mat(count(:, 2)) == 1, 1));
gene_expr = filter_rows(gene_expr, valid);

[~, order] = sort(gene_expr.rows.gene);
gene_expr = filter_rows(gene_expr, order);

fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '~/rnaseq_gene_expression');








% PROSTATE CANCER - PRENSNER (2011), PATIENTS
cd ~/datasets/prensner_pca_rnaseq/patients/reads
mkdir ../tophat_alignments
echo *.fq.gz | parallel -J prensner_align -n4 -c8 -m20 'tophat2 -o ../tophat_alignments/${x%.fq.gz} --transcriptome-index ~/tools/tophat-indexes/homo_sapiens/iGenome_37.2_NCBI ~/tools/bowtie2-indexes/homo_sapiens/hg19 $x'

cd ~/datasets/prensner_pca_rnaseq/patients/tophat_alignments
for x in */accepted_hits.bam; do mv $x ${x%/*}.bam; done
mkdir unmapped
for x in */unmapped.bam; do mv $x unmapped/${x%/*}.bam; done
for x in */prep_reads.info; do rm -r ${x%/*}; done
echo *.bam | parallel -J index_bam -n8 'samtools index $x'

sam count -b ~/organisms/homo_sapiens/ensembl_68/exons.composite.with_chr.bed *.bam







% PROSTATE CANCER - KANNAN ET AL.
cd ~/datasets/kannan_pc_rnaseq/reads
mkdir ../tophat_alignments
echo *_1.fastq.gz | parallel -Jkannan_align -n6 -c8 -m20 'tophat2 -o ../tophat_alignments/${x%_1.fastq.gz} --transcriptome-index ~/tools/tophat-indexes/homo_sapiens/iGenome_37.2_NCBI ~/tools/bowtie2-indexes/homo_sapiens/hg19 $x ${x%_1.fastq.gz}_2.fastq.gz'

cd ~/datasets/kannan_pc_rnaseq/tophat_alignments
for x in */accepted_hits.bam; do mv $x ${x%/*}.bam; done
mkdir unmapped
for x in */unmapped.bam; do mv $x unmapped/${x%/*}.bam; done
for x in */prep_reads.info; do rm -r ${x%/*}; done
echo *.bam | parallel -J index_bam -n8 'samtools index $x'

sam count -b ~/organisms/homo_sapiens/ensembl_68/exons.composite.with_chr.bed *.bam














% PROSTATE CANCER - TCGA PRAD
[data, headers] = readtable( ...
	'~/datasets/tcga_prad/rnaseq/tophat_alignments/gene_expr.txt', ...
	'Numeric', 2:10000);
gene_expr = struct;
gene_expr.meta.sample_id = headers(2:end)';
gene_expr.rows.gene = ensembl_to_hugo(data{1}, ...
	'~/organisms/homo_sapiens/ensembl_68/Homo_sapiens.GRCh37.68.gtf');
gene_expr.mean = cat(2, data{2:end}) + 1;
gene_expr.mean = normalize_median_of_ratios(gene_expr.mean, 1000);

gene_expr.meta.tcga_sample_id = gene_expr.meta.sample_id;
for k = 1:length(gene_expr.meta.sample_id)
	gene_expr.meta.sample_id{k} = gene_expr.meta.sample_id{k}(7:13);
end

save ~/datasets/tcga_prad/rnaseq/gene_expression.mat gene_expr

valid = ~rx(gene_expr.rows.gene, '^(RP11|AC\d)');
gene_expr = filter_rows(gene_expr, valid);

count = tabulate(gene_expr.rows.gene);
valid = ismember(gene_expr.rows.gene, count(cell2mat(count(:, 2)) == 1, 1));
gene_expr = filter_rows(gene_expr, valid);

[~, order] = sort(gene_expr.rows.gene);
gene_expr = filter_rows(gene_expr, order);

fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '/dev/shm/rnaseq_gene_expression');









% PROSTATE CANCER - PC CELL LINES
[data, headers] = readtable( ...
	'~/datasets/pc_celllines/rna-seq/tophat_alignments/gene_expr.txt', ...
	'Numeric', 2:10000);
gene_expr = struct;
gene_expr.meta.sample_id = headers(2:end)';
gene_expr.rows.gene = ensembl_to_hugo(data{1}, ...
	'~/organisms/homo_sapiens/ensembl_68/Homo_sapiens.GRCh37.68.gtf');
gene_expr.mean = cat(2, data{2:end}) + 1;
gene_expr.mean = normalize_median_of_ratios(gene_expr.mean, 1000);

gene_expr.meta.tcga_sample_id = gene_expr.meta.sample_id;
for k = 1:length(gene_expr.meta.sample_id)
	gene_expr.meta.sample_id{k} = gene_expr.meta.sample_id{k}(7:13);
end

save ~/datasets/tcga_prad/rnaseq/gene_expression.mat gene_expr

valid = ~rx(gene_expr.rows.gene, '^(RP11|AC\d)');
gene_expr = filter_rows(gene_expr, valid);

count = tabulate(gene_expr.rows.gene);
valid = ismember(gene_expr.rows.gene, count(cell2mat(count(:, 2)) == 1, 1));
gene_expr = filter_rows(gene_expr, valid);

[~, order] = sort(gene_expr.rows.gene);
gene_expr = filter_rows(gene_expr, order);

fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '/dev/shm/rnaseq_gene_expression');










% PROSTATE CANCER - BOVA PATIENT A21 METASTASES
load ~/datasets/bova_metastasis/rna-seq/gene_expression
valid = ~rx(gene_expr.rows.gene, '^(RP11|AC\d)');
gene_expr = filter_rows(gene_expr, valid);

count = tabulate(gene_expr.rows.gene);
valid = ismember(gene_expr.rows.gene, count(cell2mat(count(:, 2)) == 1, 1));
gene_expr = filter_rows(gene_expr, valid);

[~, order] = sort(gene_expr.rows.gene);
gene_expr = filter_rows(gene_expr, order);

fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '~/rnaseq_gene_expression');










% PROSTATE CANCER - TAMPERE (2012)
load ~/datasets/tampere_pca/rna-seq/gene_expression
gene_expr.mean = quantilenorm(gene_expr.mean + 1);
gene_expr = filter(gene_expr, [ ...
	find(rx(gene_expr.meta.sample_id, 'BPH[kt]')), ...
	find(rx(gene_expr.meta.sample_id, 'BPHn')), ...
	find(rx(gene_expr.meta.sample_id, 'PCaN')), ...
	find(rx(gene_expr.meta.sample_id, 'PCaP')), ...
	find(rx(gene_expr.meta.sample_id, 'CRPC'))]);

fmatrix = cat_fmatrix(fmatrix_gene_expression(gene_expr), ...
	fmatrix_clinical(gene_expr));
export_betastasis(fmatrix, '~/rnaseq_gene_expression');



% Ensembl gene expression
[data, headers] = readtable( ...
	'~/datasets/tampere_pca/rna-seq/tophat_alignments/gene_expr.txt', ...
	'Numeric', 2:10000);
gene_expr = struct;
gene_expr.meta.sample_id = headers(2:end)';
gene_expr.rows.gene = ensembl_to_hugo(data{1}, ...
	'~/organisms/homo_sapiens/ensembl_68/Homo_sapiens.GRCh37.68.gtf');
gene_expr.mean = cat(2, data{2:end}) + 1;
gene_expr = normalize_ratio_medians(gene_expr, 1000);

valid = ~rx(gene_expr.rows.gene, '^(RP11|AC\d)');
gene_expr = filter_rows(gene_expr, valid);

count = tabulate(gene_expr.rows.gene);
valid = ismember(gene_expr.rows.gene, count(cell2mat(count(:, 2)) == 1, 1));
gene_expr = filter_rows(gene_expr, valid);

gene_expr = filter(gene_expr, [ ...
	find(rx(gene_expr.meta.sample_id, 'BPH[kt]')), ...
	find(rx(gene_expr.meta.sample_id, 'BPHn')), ...
	find(rx(gene_expr.meta.sample_id, 'PCaN')), ...
	find(rx(gene_expr.meta.sample_id, 'PCaP')), ...
	find(rx(gene_expr.meta.sample_id, 'CRPC'))]);

fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '~/gene_expression_ensembl');



% MicroRNA expression
load ~/datasets/tampere_pca/smallrna-seq/mirna_expression_normalized
mirna_expr = filter_rows(mirna_expr, ~rx(mirna_expr.rows.name, ' '));
mirna_expr.rows.mirna = mirna_expr.rows.name;
mirna_expr.mean = mirna_expr.mean + 1;

mirna_expr.meta = tampere_pca_clinical(mirna_expr.meta);
mirna_expr = filter(mirna_expr, [ ...
	find(rx(mirna_expr.meta.sample_id, 'BPH[kt]')), ...
	find(rx(mirna_expr.meta.sample_id, 'BPHn')), ...
	find(rx(mirna_expr.meta.sample_id, 'PCaN')), ...
	find(rx(mirna_expr.meta.sample_id, 'PCaP')), ...
	find(rx(mirna_expr.meta.sample_id, 'CRPC'))]);

fmatrix = cat_fmatrix(fmatrix_mirna_expression(mirna_expr));
export_betastasis(fmatrix, '~/mirna_expression');




meta = gene_expr.meta;
primary = filter(meta, rx(meta.sample_id, 'PCa'));
export_sample_groups('~/tumor_grade_groups.json', ...
	'All tumors', filter(meta, rx(meta.sample_id, 'PCa|CRPC')), ...
	'PCa', primary, ...
	'PCa (unprogressed)', filter(meta, rx(meta.sample_id, 'PCaN')), ...
	'PCa (progressed)', filter(meta, rx(meta.sample_id, 'PCaP')), ...
	'PCa (GS > 7)', filter(primary, primary.gleason > 7), ...
	'PCa (GS = 7)', filter(primary, primary.gleason == 7), ...
	'PCa (GS < 7)', filter(primary, primary.gleason < 7), ...
	'CRPC', filter(meta, rx(meta.sample_id, 'CRPC')), ...
	'BPH', filter(meta, rx(meta.sample_id, 'BPH')));







% JAYANT RANE DATASET
cd ~/datasets/jayant_rane_mrna/raw_array_samples
raw = import_uarray_raw('Affymetrix HG U133 Plus 2.0');

load /worktmp/pipeline/platforms/affy_hg_u133_plus_2/gene_probesets
expr = uarray_expression_rma(raw, probesets);

bph_ids = {'PE512', 'PE561', 'PE574', 'PE627', 'PE662', 'PE690', 'PE693'};
bph = false(1, length(expr.meta.sample_id));
for k = 1:length(bph_ids)
	bph(rx(expr.meta.sample_id, bph_ids{k})) = true;
end

gleason = [6 7 9 7 nan nan 7 8 7 nan 6 nan 7 6 7 nan nan 7];
expr.meta.gleason_grade = nan(1, length(expr.meta.sample_id));
expr.meta.gleason_grade(1:2:end) = gleason;
expr.meta.gleason_grade(2:2:end) = gleason;

stem = rx(expr.meta.sample_id, 'CD133');
committed = rx(expr.meta.sample_id, 'Alpha2Low');

fmatrix = fmatrix_gene_expression(expr);
export_betastasis(fmatrix, ...
	'~/web/betadata/prostate_cancer/birnie08/gene_expression');

export_sample_groups('~/web/betadata/prostate_cancer/birnie08/groups.json', ...
	'Prostate cancer', filter(expr, ~bph), ...
	'Prostate cancer (GG > 6)', ...
		filter(expr, ~bph & expr.meta.gleason_grade > 6), ...
	'BPH', filter(expr, bph), ...
	'Stem cell', filter(expr, stem), ...
	'Differentiated', filter(expr, committed));










% GLIOMA - TCGA GBM
expr = realize(query('tcga/gbm affy huex gene expression'));

hms = query('gbm agilent cgh 244a hms');
mskcc = query('gbm agilent cgh 244a mskcc');
load([ppath '/platforms/agilent_hg_cgh_244a/cgh_probesets']);
hms.Sample.ID = strcat(hms.Sample.ID, ' (HMS)');
mskcc.Sample.ID = strcat(mskcc.Sample.ID, ' (MSKCC)');
combined = query_union(hms, mskcc);
tumor = filter_query(combined, 'sample type ~ tumor');
normal = filter_query(combined, 'sample type ~ normal|reference');
[tumor, normal] = paired_samples(tumor, normal, 'Filename')
tumor = realize(tumor);
normal = realize(normal);
[~, order] = sort(tumor.Meta.Sample.ID);
tumor = filter_query(tumor, order);
normal = filter_query(normal, order);
gene_lr = gene_cgh_logratios(tumor, normal, probesets, 'MinRegion', 100e3);

% Kaplan-Meier plots for protein expression data
load ~/datasets/tcga_gbm/protein_expression
fmatrix = cat_fmatrix(fmatrix_protein_expression(prot_expr), ...
	fmatrix_clinical(prot_expr));

export_betastasis(fmatrix, '~/protein_expression');
export_sample_groups('~/protein_expression_groups.json', ...
	'All tumors', prot_expr);














% COLORECTAL CANCER - TCGA
exon_expr = realize(query('tcga/crc wt exon expression new'));
exon_expr.Mean = quantilenorm(exon_expr.Mean);

fmatrix = fmatrix_exon_expression(exon_expr);
export_betastasis(fmatrix, '~/web/betadata/colorectal_cancer/tcga_crc_splice_new');

qset = query('tcga/crc wt exon expression new');
coad = filter_query(qset, 'sample type ~ colon.*adenocarcinoma');
read = filter_query(qset, 'sample type ~ rectal.*adenocarcinoma');

metastatic = filter_query(qset, strcmp(qset.Misc.distant_metastasis_pathologic_spread, 'M1'));
non_metastatic = filter_query(qset, strcmp(qset.Misc.distant_metastasis_pathologic_spread, 'M0'));

stage_i = filter_query(qset, rx(qset.Misc.tumor_stage, 'Stage I[^I]*$'));
stage_ii = filter_query(qset, rx(qset.Misc.tumor_stage, 'Stage II[^I]*$'));
stage_iii = filter_query(qset, rx(qset.Misc.tumor_stage, 'Stage III[^I]*$'));

export_sample_groups( ...
	'~/web/betadata/colorectal_cancer/tcga_crc_groups/all_groups.json', ...
	'All tumors', qset, ...
	'Colon adenocarcinomas', coad, ...
	'Rectal adenocarcinomas', read, ...
	'Metastasized', metastatic, ...
	'Non-metastasized', non_metastatic, ...
	'Stage I', stage_i, ...
	'Stage II', stage_ii, ...
	'Stage III', stage_iii);


% Analyze differential expression
exon_expr = realize(query('tcga/crc wt exon expression new'));
exon_expr.Mean = quantilenorm(exon_expr.Mean);

metastatic = filter_query(exon_expr, ...
	strcmp(qset.Misc.distant_metastasis_pathologic_spread, 'M1'));
non_metastatic = filter_query(exon_expr, ...
	strcmp(qset.Misc.distant_metastasis_pathologic_spread, 'M0'));
	
find_diff_splicing_rnaseq(metastatic, non_metastatic);


% Test differential expression detection using exon microarrays
exon_expr = realize(query('taylor prostate/huex exon expression'));
exon_expr.Mean = quantilenorm(exon_expr.Mean);

fmatrix = fmatrix_exon_splice_or(exon_expr);
export_betastasis(fmatrix, '~/web/betadata/prostate_cancer/taylor_splice_new');

find_diff_splicing_rnaseq(exon_expr, exon_expr);










% PROSTATE CANCER - FRIEDLANDER ET AL.
cd ~/datasets/friedlander_crpc_methylation/acgh_raw
raw = import_uarray_raw('Agilent HG CGH 244A');
tumor = filter(raw, rx(raw.meta.uarray_channel, 'Cy3'));  % Liver metastasis
normal = filter(raw, rx(raw.meta.uarray_channel, 'Cy5'));
[tumor, normal] = pair(tumor, normal, 'uarray_filename');

load /worktmp/pipeline/platforms/agilent_hg_cgh_244a/cgh_probesets
segments = cgh_segment_rcbs(tumor, normal, probesets)

cn_seg_to_track(segments, '~/datasets/friedlander_crpc_methylation/acgh_raw/tumor_vs_ref.seg');






% Applied Biosystems microRNA expression dataset (40 tissues)
cd ~/datasets/abi_mirna_tissue_expr
[data, headers] = readtable('microrna_in_tissues.txt');

expr = struct;
expr.rows.mirna = data{1};
expr.meta.sample_id = headers(2:end)';
expr.mean = nan(length(expr.rows.mirna), length(expr.meta.sample_id));

for s = 1:length(expr.meta.sample_id)
	expr.mean(:, s) = str2double(data{s+1});
end

expr.mean = 2.^(35 - expr.mean);

fmatrix = fmatrix_mirna_expression(expr);
export_betastasis(fmatrix, '~/abi_mirna_tissue_expr');








% PANCREATIC CANCER - CANCER CELL LINE ENCYCLOPEDIA
load ~/datasets/betastasis/cancer_cell_line_encyclopedia/gene_expression
gene_expr.meta.tissue = regexprep(gene_expr.meta.sample_id, '^.*?_', '');
gene_expr.meta.sample_id = regexprep(gene_expr.meta.sample_id, '_.*', '');
gene_expr = filter(gene_expr, ~strcmp('TT', gene_expr.meta.sample_id));

for k = 1:length(gene_expr.meta.tissue)
	gene_expr.meta.tissue{k}(2:end) = lower(gene_expr.meta.tissue{k}(2:end));
	gene_expr.meta.tissue{k} = strrep(gene_expr.meta.tissue{k}, '_', ' ');
end
tissue_types = unique(gene_expr.meta.tissue);

gene_expr.mean = 2.^gene_expr.mean;
fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '~/ccle_gene_expression');

group_params = {};
meta = gene_expr.meta;
for k = 1:length(tissue_types)
	group_params{end+1} = tissue_types{k};
	group_params{end+1} = filter(meta, strcmp(tissue_types{k}, meta.tissue));
end

export_sample_groups('~/tissue_groups.json', group_params{:});










% PANCREATIC CANCER - SHAIN ET AL. (2011)
load ~/datasets/betastasis_datasets/rauhala_betastasis/gene_expression
data = readtable('/data/csb/platforms/agilent_hg_4x44k/GPL4133-9097.txt', ...
	'Comment', '^#');
featid_to_gene = containers.Map(data{1}, data{10});
gene_expr.rows.gene_symbol = featid_to_gene.values(gene_expr.rows.id);
valid = ~cellfun(@isempty, gene_expr.rows.gene_symbol);
gene_expr = filter_rows(gene_expr, valid);
featid_to_tx = containers.Map(data{1}, data{8});
gene_expr.rows.transcript = featid_to_tx.values(gene_expr.rows.id);
valid = ~cellfun(@isempty, gene_expr.rows.transcript);
gene_expr = filter_rows(gene_expr, valid);

[uniq_tx, uidx] = unique(gene_expr.rows.transcript);
ps_mean = nan(length(uniq_tx), size(gene_expr.mean, 2));
for k = 1:length(uniq_tx)
	ps_mean(k, :) = mean(gene_expr.mean( ...
		strcmp(gene_expr.rows.transcript, uniq_tx{k}), :), 1);
end

gene_expr = filter_rows(gene_expr, uidx);
gene_expr.mean = ps_mean;
gene_expr.mean = 2.^gene_expr.mean;

xenografts = find(rx(gene_expr.meta.sample_source_ch1, 'xenograft'));
cell_lines = find(rx(gene_expr.meta.sample_source_ch1, 'cell line'));
gene_expr = filter(gene_expr, [ xenografts, cell_lines ]);

fmatrix = fmatrix_transcript_expression(gene_expr);
export_betastasis(fmatrix, '/data/tmp/annalam/shain_et_al_2011');
[~, ~] = mkdir('/data/tmp/annalam/shain_et_al_2011/groups');
export_sample_groups( ...
	'/data/tmp/annalam/shain_et_al_2011/groups/cell_xenograft_groups.json', ...
	'All samples', gene_expr, ...
	'Xenografts', filter(gene_expr, xenografts), ...
	'Cell lines', filter(gene_expr, cell_lines));


cd ~/datasets/betastasis_datasets/rauhala_betastasis/acgh
raw = import_uarray_raw('Agilent 244A');
load /data/csb/pipeline/platforms/agilent_hg_cgh_244a/cgh_probesets
geo = read_geo_series_matrix('../GSE25273_series_matrix.txt');
geo.meta.sample_id = geo.meta.sample_title;
tumor = filter(raw, rx(raw.meta.uarray_channel, 'Cy5'));
normal = filter(raw, rx(raw.meta.uarray_channel, 'Cy3'));

[~, pos] = ismember(strrep(tumor.meta.sample_id, '(Cy5)', ''), geo.meta.sample_accession);
tumor.meta.sample_id = geo.meta.sample_title(pos);
tumor.meta.sample_type = geo.meta.sample_source_ch1(pos);

segments = cgh_segment_rcbs(tumor, normal, probesets);

xenografts = find(rx(tumor.meta.sample_type, 'xenograft'));
cell_lines = find(rx(tumor.meta.sample_type, 'cell line'));
segments = filter(segments, [ xenografts, cell_lines ]);

cn_seg_to_track(segments, '~/shain_et_al_2011_acgh_cbs.seg');
!scp ~/shain_et_al_2011_acgh_cbs.seg intianjora.cs.tut.fi:~/betadata/pancreatic_cancer/shain_et_al_2011/acgh_cbs.seg







% MAURO MICROARRAY DATASET
cd ~/datasets/betastasis/mauro_bladder_cancer
load raw
S = size(raw.mean, 2);
raw.meta.sample_id = repmat({'Reference RNA'}, 1, S);
raw.meta.sample_id(2:2:S) = { ...
	'SCaBER', '5637', ...
	'T24 #1', 'J82', 'HT-1376', 'UM-UC-3', ...
	'T24 #2', 'TCCSUP #1', '143 #1', '63', ...
	'TCCSUP #2', '143 #2', '88', '185' };

load /data/csb/pipeline/platforms/agilent_human_ge_v1/probesets
gene_expr = uarray_expression_rma(raw, probesets);
save gene_expression.mat gene_expr

load gene_expression
gene_expr.rows.gene = organism.Genes.Name;
gene_expr = filter(gene_expr, 2:2:S);
%for s = 1:size(gene_expr.mean, 2)/2
%	diff_expr.mean(:, s) = ...
%		log2(gene_expr.mean(:, s*2) ./ gene_expr.mean(:, s*2-1));
%end

[~, order] = sort(gene_expr.meta.sample_id);
gene_expr = filter(gene_expr, order);

fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '~/mauro_expression');
!mkdir ~/mauro_expression/groups
export_sample_groups('~/mauro_expression/groups/groups.json', ...
	'All samples', gene_expr);









% EXPO DATASET
load ~/datasets/intgen_expo/gene_expression
gene_expr = filter(gene_expr, rx(gene_expr.meta.organ, 'prostate'));
fmatrix = fmatrix_gene_expression(gene_expr);
export_betastasis(fmatrix, '~/expo_pca');
export_sample_groups('~/groups.json', ...
	'All samples', gene_expr);
!mkdir ~/expo_pca/groups && mv ~/groups.json ~/expo_pca/groups/





% PROSTATE CANCER - TAMURA (2007)
cd ~/datasets/betastasis_datasets/tamura07
geo = read_geo_series_matrix('GSE6811_series_matrix.txt', 'ReadData', true);
[data, headers] = readtable('GPL4747.annot', 'Comment', '^[\^!#]');
for k = 1:length(data), data{k} = data{k}(1:end-1); end

probe_id_to_gene = containers.Map(data{1}, data{3});
geo.rows.gene_symbol = probe_id_to_gene.values(data{1});
valid = ~cellfun(@isempty, geo.rows.gene_symbol);

fmatrix = fmatrix_gene_expression(diff_expr);
export_betastasis(fmatrix, '~/tamura07');
!mkdir ~/tamura07/groups
export_sample_groups('~/tamura07/groups/groups.json', ...
	'All samples', diff_expr);





% GLIOMA - PARKER ET AL (2012)
load ~/datasets/wei_glioma_seq/gene_expression
norm_expr = normalize_rnaseq_expr(gene_expr, 'RPKM');
fmatrix = fmatrix_gene_expression(norm_expr);
export_betastasis(fmatrix, '~/wei_glioma_seq');

