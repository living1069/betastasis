%% datadir  '/example/testdir'
%% platform 'agilent_244k_tcga_custom_03','affy_huex_1_0_st','affy_ht_hg_u133a'
%% datasetname 'testA' 
%% pipelinedir '/example/pipelinedir'

function smrTest(datadir,pipelinedir,platform,datasetname)
%% this is simply from the 'start' script
currentdir=pwd
cd ([pipelinedir '/sources']); 
pipeline;

%% next CD to the data directory for the HuEx data ...
cd (datadir)

platform_path = [pipelinedir '/platforms/' platform];

if (strcmp(platform,'agilent_244k_tcga_custom_03'))
	platformname='Agilent 244K TCGA custom 3';
    	dge_affinities_str = 'agi3_to_dge_affinities';
	dge_affinities_path = [platform_path '/' dge_affinities_str '.05sep10.mat'];
elseif (strcmp(platform,'affy_huex_1_0_st'))
	platformname='Affymetrix Human Exon 1.0 ST';
    	dge_affinities_str = 'huex_to_dge_affinities';
	dge_affinities_path = [platform_path '/' dge_affinities_str '.mat'];
elseif (strcmp(platform,'affy_ht_hg_u133a'))
	platformname='Affymetrix HT HG U133A';
    	dge_affinities_str = 'ht_hg_u133a_to_dge_affinities';
	dge_affinities_path = [platform_path '/' dge_affinities_str '.mat'];
end

remove_dataset(datasetname);
import_tcga_uarray_data ( datasetname, platformname )

data = realize ( query ( datasetname ) )

%% at this point, data is a struct with
%%	data.Mean	223302 x 1052 double
%%	data.Meta.Type	'Microarray probe intensities'
%%	data.Meta.Sample.ID		1052x1 cell	
%%	data.Meta.Sample.Channel	1052x1 cell
%%	data.Meta.Sample.Filename	1052x1 cell
%%
%% note, even #s appear to have TCGA barcodes, and odd #s have the 'reference'
%% also, the same filename is given for both channels

gene_probeset_path = [platform_path '/gene_probesets.mat'];

gene_probeset_name_str = [platform '_gene_probesets'];

load (gene_probeset_path)
gene_probeset_name_var = eval([gene_probeset_name_str]);

expr = uarray_expression_rma ( data, gene_probeset_name_var )

load (dge_affinities_path)
dge_affinities_var = eval([dge_affinities_str]);
expr.Mean = expr.Mean ./ repmat ( dge_affinities_var, 1, size(expr.Mean, 2) )

resultsfile = [strrep(datasetname,' ','_') '_new_expr.mat'];
cd (currentdir)
save (resultsfile,'expr')
