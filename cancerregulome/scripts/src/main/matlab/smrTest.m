%% datadirectory '/titan/cancerregulome3/TCGA/clinical-data-repository/tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distroymous/tumor/ov/cgcc/unc.edu/agilentg4502a_07_3/transcriptome/'
%% platform 'agilent_244k_tcga_custom_03','affy_huex_1_0_st','affy_ht_hg_u133a'
%% datasetname 'g4502a_07_3 testA' 
function smrTest(datadirectory,platform,datasetname)
%% this is simply from the 'start' script
cd /titan/cancerregulome5/csbPipeline/sources; 
pipeline;

%% next CD to the data directory for the HuEx data ...
cd (datadirectory)

platform_path = ['/titan/cancerregulome5/csbPipeline/platforms/' platform];

if (strcmp(platform,'agilent_244k_tcga_custom_03'))
	platformname='Agilent 244K TCGA custom 3';
    dge_affinities = 'agi3_to_dge_affinities';
	dge_affinities_path = [platform_path '/' dge_affinities '.05sep10.mat'];
elseif (strcmp(platform,'affy_huex_1_0_st'))
	platformname='Affymetrix Human Exon 1.0 ST';
    dge_affinities = 'huex_to_dge_affinities';
	dge_affinities_path = [platform_path '/' dge_affinities '.mat'];
elseif (strcmp(platform,'affy_ht_hg_u133a'))
	platformname='Affymetrix HT HG U133A';
    dge_affinities = 'ht_hg_u133a_to_dge_affinities';
	dge_affinities_path = [platform_path '/' dge_affinities '.mat'];
end


import_tcga_uarray_data ( datasetname, platformname )

data = realize ( query ( datasetname ) )

%% at this point, data is a struct with
%%	data.Mean	223302 x 1052 double
%%	data.Meta.Type	'Microarray probe intensities'
%%	data.Meta.Sample.ID		1052x1 cell	eg TCGA-29-1697-01 or Stratagene Universal Reference
%%	data.Meta.Sample.Channel	1052x1 cell	eg Cy5		   or Cy3
%%	data.Meta.Sample.Filename	1052x1 cell	eg unc.edu_OV.AgilentG4502A_07_3.Level_1.1.5.0/US82800149_251976010696_S01_GE2-v5_10_Apr08.txt
%%
%% note, even #s appear to have TCGA barcodes, and odd #s have the 'reference'
%% also, the same filename is given for both channels

gene_probeset_path = [platform_path '/gene_probesets.mat'];

gene_probeset_name = [platform '_gene_probesets'];

load gene_probeset_path

expr = uarray_expression_rma ( data, gene_probeset_name )

load dge_affinities_path
expr.Mean = expr.Mean ./ repmat ( dge_affinities, 1, size(expr.Mean, 2) )

resultsfile = [strrep(datasetname,' ','_') '_new_expr.mat'];
save resultsfile expr