from __future__ import print_function
import os, sys, json

import bottle
from collections import OrderedDict
from bottle import route, template, static_file, abort, redirect

hierarchy = OrderedDict()

def name_to_url(name):
	name = name.lower()
	name = name.replace('(', '')
	name = name.replace(')', '')
	name = name.replace(' ', '_')
	name = name.replace('.', '')
	name = name.replace('#', '')
	name = name.replace('_>_', '/')  # FIXME: Ignore amount of whitespace
	return name

def btemplate(tmpl, name, **kwargs):
	kwargs.setdefault('category', '')
	hierarchy[name_to_url(name)] = (name, tmpl, kwargs)

def disease(name, **kwargs): btemplate('disease', name, **kwargs)
def dataset(name, **kwargs): btemplate('dataset', name, **kwargs)
def gene_boxplot(name, **kwargs): btemplate('gene_boxplot', name, **kwargs)
def kaplan_meier(name, **kwargs): btemplate('kaplan_meier', name, **kwargs)
def expr_barplot(name, **kwargs): btemplate('expr_barplot', name, **kwargs)
def exon_expr_barplot(name, **kwargs):
	btemplate('exon_expr_barplot', name, **kwargs)
def feature_matrix(name, **kwargs): btemplate('feature_matrix', name, **kwargs)
def gene_scatter(name, **kwargs): btemplate('gene_scatter', name, **kwargs)
def aberrated_scatter(name, **kwargs):
	btemplate('aberrated_scatter', name, **kwargs)
def parallel_coordinates(name, **kwargs):
	btemplate('parallel_coordinates', name, **kwargs)
def geneset_correlation(name, **kwargs):
	btemplate('gset_correlation', name, **kwargs)
def alt_splicing(name, **kwargs): btemplate('splice_boxplot', name, **kwargs)
def alt_splicing_exons(name, **kwargs):
	btemplate('splice_exon_boxplot', name, **kwargs)
def alt_splicing_survival(name, **kwargs):
	btemplate('splice_survival', name, **kwargs)
def timeseries(name, **kwargs): btemplate('timeseries', name, **kwargs)
def graph(name, **kwargs): btemplate('graph', name, **kwargs)
def mutation_list(name, **kwargs): btemplate('mutation_list', name, **kwargs)
def cdna_uarray_probes(name, **kwargs):
	btemplate('cdna_uarray_probes', name, **kwargs)
def igv_track(name, **kwargs): btemplate('igv_track', name, **kwargs)
def circvis(name, **kwargs): btemplate('circvis', name, **kwargs)

@route('/js/:filename#.+#')
def get_script(filename):
	return static_file(filename, root='js')

@route('/style/:filename#.+#')
def get_style(filename):
	return static_file(filename, root='style')

@route('/images/:filename#.+#')
def get_image(filename):
	return static_file(filename, root='images')

@route('/data/:filename#.+#')
def redirect_data(filename):
	return static_file(filename, root='/home/csbgroup/public_htmlhttp://betastasis.cs.tut.fi/betadata')

@route('/')
def index():
	return template('index', hierarchy=hierarchy, path='', url='')

@route('/:page#.+#')
def get_page(page):
	if page[-1] == '/': page = page[:-1]
	if not page in hierarchy:
		abort(404, 'Page %s does not exist.' % page)
	
	path, tmpl, kwargs = hierarchy[page]
	kwargs['hierarchy'] = hierarchy
	kwargs['path'] = path
	kwargs['url'] = page
	return template(tmpl, **kwargs)











disease('Bladder cancer')

disease('Breast cancer')

disease('Colorectal cancer')

#disease('Endometrial cancer')

disease('Glioma', content='''
<p>A glioma is a type of tumor that arises from glial cells. Gliomas usually occur in the brain, but sometimes also in the spine.</p>''')

disease('Pancreatic cancer')

disease('Prostate cancer', content='''
<p>Prostate cancer is a form of cancer that develops in the prostate, a gland in the male reproductive system. Most prostate cancers are slow growing; however, there are cases of aggressive prostate cancers. The cancer cells may metastasize from the prostate to other parts of the body, particularly the bones and lymph nodes. Prostate cancer usually progresses from adenocarcinoma in situ to adenocarcinoma and then metastatic disease.</p>''')

disease('Tissues', contents='''
<p>This section contains datasets that characterize RNA expression, protein levels or epigenetic features across multiple tissue types. Only normal tissues are included in these datasets.</p>
''')

disease('Platforms', contents='''
<p>This section contains technical information about microarray and high throughput sequencing platforms.</p>
''')

disease('Tools', content='''
<p>This section contains tools to upload and research your own data.</p>
''')









##################
# BLADDER CANCER #
##################


# Scaravilli (2012)
dataset('Bladder cancer > Scaravilli (2012)',
	title='Scaravilli - Bladder cancer cell lines and patients',
	stats='', hidden=True,
	platform='Agilent Whole Human Genome 4x44K')

expr_barplot('Bladder cancer > Scaravilli (2012) > Gene expression barplot',
	category='Gene expression microarrays', gene='TP53',
	data='http://betastasis.cs.tut.fi/betadata/bladder_cancer/scaravilli_2012/gene_expression',
	platform='Agilent Whole Human Genome 4x44K')



# TCGA - Bladder cancer (BLCA)
dataset('Bladder cancer > TCGA BLCA',
	title='TCGA - Bladder urothelial carcinoma (BLCA)',
	stats='Total samples: 67<br>',
	platform='Illumina HiSeq 2000')

expr_barplot('Bladder cancer > TCGA BLCA > Gene expression barplot',
	category='Whole transcriptome sequencing', gene='TP53',
	data='http://betastasis.cs.tut.fi/betadata/bladder_cancer/tcga_prad',
	platform='Illumina HiSeq 2000')
















#################
# BREAST CANCER #
#################

# TCGA - Breast cancer (BRCA)
dataset('Breast cancer > TCGA BRCA',
	title='TCGA - Breast cancer (BRCA)',
	stats='Total samples: 471<br>Total patients: 435<br>',
	platform='<a href="/platforms/cdna_microarrays/agilent_244k_tcga_custom_3">Agilent 244K TCGA custom #3</a><br>Affymetrix Genome Wide SNP 6.0<br>')

gene_boxplot('Breast cancer > TCGA BRCA > Gene expression box plot',
	category='Gene-level views', gene='TP53',
	data='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_expr',
	groups='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_groups/t_stage_boxplot.json',
	platform='<a href="/platforms/cdna_microarrays/agilent_244k_tcga_custom_3">Agilent 244K TCGA custom 3</a>')

feature_matrix('Breast cancer > TCGA BRCA > Feature matrix',
	category='Gene-level views', data_root='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_cna', groups='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_groups/t_stage.json', platform='Affymetrix GW SNP 6.0', sample_labelsize='70', default_features='TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),ERBB2 (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)')

kaplan_meier('Breast cancer > TCGA BRCA > Kaplan-Meier survival plot',
	category='Gene-level views',
	data='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_expr',
	groups='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_groups/t_stage.json',
	feature='TP53',
	platform='<a href="/platform/cdna_microarrays/agilent_244k_tcga_custom_3">Agilent 244K TCGA custom 3</a>')

parallel_coordinates('Breast cancer > TCGA BRCA > Parallel coordinates',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_expr',
	groups='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_groups/t_stage.json',
	platform='<a href="/platforms/cdna_microarrays/agilent_244k_tcga_custom_3">Agilent 244K TCGA custom 3</a>',
	features='["ERBB2", "PTEN", "TP53", "PDGFRA", "BRCA1"]')

igv_track('Breast cancer > TCGA BRCA > All tumors, paired with adjacent normal (Affymetrix GW SNP 6.0)',
	category='Copy number tracks',
	data='http://betastasis.cs.tut.fi/betadata/breast_cancer/tcga_tracks/tumors_vs_adj_normal.seg')
















#####################
# COLORECTAL CANCER #
#####################

# Veeriah et al. (2009)
dataset('Colorectal cancer > Veeriah et al. (2009)',
	title='Veeriah et al. - Somatic mutations of the Parkinson''s disease-associated gene PARK2 in glioblastoma and other human malignancies',
	pubmed_id='19946270',
	publication_date='29-11-2009',
	stats='Total samples: 98<br>Total patients: 95<br>',
	platform='Agilent HG CGH 244A<br>',
	content='''
<p>Agilent 244K aCGH data for 95 human colon cancer samples and 3 colon cancer cell lines.</p>
''')

feature_matrix('Colorectal cancer > Veeriah et al. (2009) > Feature matrix',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/colorectal_cancer/taylor_cna',
	groups='http://betastasis.cs.tut.fi/betadata/colorectal_cancer/taylor_groups/groups.json',
	platform='Agilent HG CGH 244A', sample_labelsize='70',
	default_features='TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)')

igv_track('Colorectal cancer > Veeriah et al. (2009) > All tumors, paired with universal reference DNA',
	category='Copy number tracks', data='http://betastasis.cs.tut.fi/betadata/crc_taylor_tracks/taylor_crc.seg')
	


	

# TCGA - Colorectal cancer
dataset('Colorectal cancer > TCGA COAD & READ',
	stats='Last update: 29-06-2011<br>Total samples: 404<br>Total patients: 202<br>',
	platform='Illumina Genome Analyzer II<br>Affymetrix Genome Wide SNP 6.0<br>',
	content='''
<p>In this dataset, we have combined the TCGA colon adenocarcinoma (COAD) dataset with the TCGA rectal adenocarcinoma (READ) dataset.</p>
''')

igv_track('Colorectal cancer > TCGA COAD & READ > All COAD tumors, paired with Promega universal reference DNA',
	category='Copy number tracks',
	data='http://betastasis.cs.tut.fi/betadata/colorectal_cancer/tcga_coad_tracks/tumor_vs_ref.seg.gz')











##########
# GLIOMA #
##########

dataset('Glioma > Parker et al. (2012)',
	title='Parker et al. (2012) - FGFR3-TACC3 fusion escapes miR-99a regulation and promotes tumorigenesis in glioblastoma',
	hidden=True,
	platform='ABI SOLiD 3')

expr_barplot('Glioma > Parker et al. (2012) > Gene expression barplot',
	category='Transcriptomic analysis', gene='TP53',
	data='http://betastasis.cs.tut.fi/betadata/glioma/parker_et_al_2012',
	platform='ABI SOLiD 3', description='''
<p>10 sample pools were sequenced using the whole transcriptome sequencing protocol of the ABI SOLiD 3 platform. RPKM expression values were calculated from raw read counts aligning to each gene.</p>
''')





# REMBRANDT
dataset('Glioma > REMBRANDT',
	title='REMBRANDT - REpository for Molecular BRAin Neoplasia DaTa',
	stats='Total samples: 524<br>Total patients: 524<br>',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 v2.0 Plus</a><br>',
	content='''
<p>REpository for Molecular BRAin Neoplasia DaTa (REMBRANDT) is a robust bioinformatics knowledgebase framework that leverages data warehousing technology to host and integrate clinical and functional genomics data from clinical trials involving patients suffering from gliomas.</p>
''')

gene_boxplot('Glioma > REMBRANDT > Gene expression in glioma subtypes',
	category='Gene-level views', gene='TP53', data='http://betastasis.cs.tut.fi/betadata/rembrandt_foo',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 v2.0 Plus</a>')

kaplan_meier('Glioma > REMBRANDT > Kaplan-Meier survival curve',
	category='Gene-level views',
	data='http://betastasis.cs.tut.fi/betadata/rembrandt_foo',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/rembrandt_groups/all_groups.json',
	feature='EGFR',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 v2.0 Plus</a>')

gene_scatter('Glioma > REMBRANDT > Two-gene scatterplot',
	category='Gene-level views', gene_x='EGFR', gene_y='GRB2',
	data='http://betastasis.cs.tut.fi/betadata/rembrandt_foo',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/rembrandt_groups/all_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_ht_hg_u133a/">Affymetrix HT HG U133A</a>')

parallel_coordinates('Glioma > REMBRANDT > Parallel coordinates',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/rembrandt_foo',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/rembrandt_groups/all_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 v2.0 Plus</a>',
	features='["EGFR", "PTEN", "TP53", "PDGFRA", "IDH1"]')

geneset_correlation('Glioma > REMBRANDT > Gene set correlation',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/rembrandt_foo',
	scatter_url='/glioma/rembrandt/gene_scatter',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/rembrandt_groups/all_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 v2.0 Plus</a>',
	ref_gene='IGFBP2', target_genes='MYC,ITGA5,PTEN,TP53,EGFR,STAT3')






# TCGA - Glioblastoma
dataset('Glioma > TCGA GBM',
	title='TCGA - Glioblastoma',
	stats='Total samples: 454<br>',
	platform='''
<a href="/platforms/cdna_microarrays/affymetrix_ht_hg_u133a/">Affymetrix HT HG U133A</a><br>
<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a><br>
Affymetrix Human miRNA 8x15K v1<br>
Affymetrix Human miRNA 8x15K v2<br>
Agilent HG CGH 244A<br>''',
	content='''
<p>The TCGA GBM dataset consists of primary tumor samples from roughly 400 patients. The dataset contains data from gene expression, CGH and methylation microarrays.</p>
''')

gene_boxplot('Glioma > TCGA GBM > Expression box plot (Affymetrix HT HG U133A)',
	category='Transcriptome analysis', gene='TP53', data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_u133',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/boxplot_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_ht_hg_u133a/">Affymetrix HT HG U133A</a>')

gene_boxplot(
	'Glioma > TCGA GBM > Expression box plot (Affymetrix Human Exon 1.0 ST)',
	category='Transcriptome analysis', gene='TP53', data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_huex',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/boxplot_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affy Human Exon 1.0 ST</a>')

gene_boxplot(
	'Glioma > TCGA GBM > Expression box plot (Agilent Human miRNA 8x15K)',
	category='Transcriptome analysis', gene='hsa-miR-21',
	data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_mirna_expr',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/boxplot_groups.json',
	platform='Agilent Human miRNA 8x15K')

gene_scatter(
	'Glioma > TCGA GBM > Two-gene scatterplot (Affymetrix HT HG U133A)',
	category='Transcriptome analysis', gene_x='EGFR', gene_y='GRB2',
	data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_u133',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_ht_hg_u133a/">Affymetrix HT HG U133A</a>')

gene_scatter(
	'Glioma > TCGA GBM > Two-gene scatterplot (Affymetrix Human Exon 1.0 ST)',
	category='Transcriptome analysis', gene_x='EGFR', gene_y='GRB2',
	data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_huex',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

geneset_correlation('Glioma > TCGA GBM > Gene set correlation',
	category='Transcriptome analysis', data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_u133',
	scatter_url='/glioma/tcga_gbm/gene_scatter_u133',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 v2.0 Plus</a>',
	ref_gene='EGFR', target_genes='MYC,ITGA5,PTEN,TP53,STAT3')

feature_matrix('Glioma > TCGA GBM > Feature matrix',
	category='Integrative analysis', data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_cna',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	platform='Agilent HG CGH 244A', sample_labelsize='100',
	default_features='EGFR (CNA),PDGFRA (CNA),TP53 (CNA),PTEN (CNA),MIIP (CNA),IGFBP2 (CNA),MTOR (CNA),AKT1 (CNA)')

parallel_coordinates('Glioma > TCGA GBM > Parallel coordinates',
	category='Integrative analysis', data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_u133',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_ht_hg_u133a/">Affymetrix HT HG U133A</a><br>Agilent Human miRNA 8x15K v1<br>Agilent Human miRNA 8x15K v2',
	features='["EGFR", "PTEN", "TP53", "PDGFRA", "IDH1"]')

kaplan_meier('Glioma > TCGA GBM > Gene survival association (Affymetrix HT HG U133A',
	category='Survival analysis',
	data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_u133',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	feature='EGFR',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_ht_hg_u133a/">Affymetrix HT HG U133A</a>')

kaplan_meier('Glioma > TCGA GBM > Gene survival association (Affymetrix Human Exon 1.0 ST)',
	category='Survival analysis',
	data='http://betastasis.cs.tut.fi/betadata/tcga_gbm_huex',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	feature='EGFR',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affy Human Exon 1.0 ST</a>')

kaplan_meier(
	'Glioma > TCGA GBM > MicroRNA survival association (Agilent Human miRNA 8x15K)',
	category='Survival analysis', data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_mirna_expr',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	feature='hsa-miR-21', feature_type='miRNA',
	platform='Agilent Human microRNA 8x15K')

kaplan_meier('Glioma > TCGA GBM > Protein survival association',
	category='Survival analysis',
	data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm/protein_expression',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm/groups/protein_expression.json',
	platform='RPPA', feature='Akt_pS473-R-V', feature_type='Protein')

alt_splicing(
	'Glioma > TCGA GBM > Alternative splicing boxplot (gene)',
	category='Alternative splicing', gene='PDGFRA',
	groups='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_groups/all_groups.json',
	data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_splice',
	default_test_group='All tumors', default_ref_group='Normal',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

alt_splicing_exons(
	'Glioma > TCGA GBM > Alternative splicing boxplot (exon)', 
	category='Alternative splicing', gene='PDGFRA',
	data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_splice',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

alt_splicing_survival(
	'Glioma > TCGA GBM > Alternative splicing survival association',
	category='Alternative splicing', gene='PDGFRA',
	data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_splice',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

igv_track('Glioma > TCGA GBM > All tumors, paired with Promega universal reference DNA',
	category='Copy number tracks', data='http://betastasis.cs.tut.fi/betadata/glioma/tcga_gbm_tracks/agilent_244a_tumor_vs_ref.seg.gz')








#####################
# PANCREATIC CANCER #
#####################

dataset('Pancreatic cancer > Shain et al. (2011)',
	title='Convergent structural alterations define SWI/SNF chromatin remodeler as a central tumor suppressive complex in pancreatic cancer',
	stats='Total samples: 68<br>',
	platforms='Agilent HG 4x44K<br>Agilent HG CGH 244A<br>',
	content='''
<p>Expression profiling and array comparative genomic hybridization of pancreatic xenografts and cell lines.</p>
''')

expr_barplot('Pancreatic cancer > Shain et al. (2011) > Gene expression barplot',
	category='Transcriptomic analysis', gene='TP53 (NM_000546)',
	data='http://betastasis.cs.tut.fi/betadata/pancreatic_cancer/shain_et_al_2011',
	platform='Agilent HG 4x44K', description='''
<p>Barplot of the expression of a single transcript in multiple samples. Pooled RNA from 11 varied cancer cell lines was used as control. Differential gene expression values imported directly from the GEO series matrix. Samples from JH21 to 235 are xenografts, while the samples from HPDE to Panc05.04 are cell lines.</p>
''')

gene_scatter('Pancreatic cancer > Shain et al. (2011) > Two-gene scatterplot',
	category='Transcriptomic analysis', data='http://betastasis.cs.tut.fi/betadata/pancreatic_cancer/shain_et_al_2011',
	gene_x='CDKN2A (NM_058197)', gene_y='CEACAM6 (NM_002483)',
	groups='http://betastasis.cs.tut.fi/betadata/pancreatic_cancer/shain_et_al_2011/groups/cell_xenograft_groups.json',
	platform='Agilent HG 4x44K')

#geneset_correlation('Pancreatic cancer > Shain et al. (2011) > Gene set correlation',
#	category='Transcriptomic analysis', data='http://betastasis.cs.tut.fi/betadata/pancreatic_cancer/shain_et_al_2011',
#	scatter_url='/pancreatic_cancer/shain_et_al_2011/two-gene_scatterplot',
#	groups='http://betastasis.cs.tut.fi/betadata/pancreatic_cancer/shain_et_al_2011/groups/cell_xenograft_groups.json',
#	platform='Agilent HG 4x44K', ref_gene='CDKN2A (NM_058197)',
#	target_genes='CEACAM6 (NM_002483)')

igv_track('Pancreatic cancer > Shain et al. (2011) > All tumors, paired with pooled normal genomic DNA from 8 individuals',
	category='Copy number tracks', data='http://betastasis.cs.tut.fi/betadata/pancreatic_cancer/shain_et_al_2011/acgh_cbs.seg')
















###################
# PROSTATE CANCER #
###################


# Prostate cancer - Bova patient A21
dataset('Prostate cancer > Bova patient A21',
	title='Bova patient A21', hidden=True,
	stats='',
	platforms='Illumina HiSeq 2000<br>',
	content='''
<p>High throughput sequencing of prostate cancer metastases from patient A21. Whole transcriptome sequencing, whole genome sequencing, and methylated DNA immunoprecipitation sequencing data is included.</p>
''')

expr_barplot('Prostate cancer > Bova patient A21 > Gene expression barplot',
	category='Whole transcriptome analysis', gene='AR',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/bova_patient_a21/rnaseq_gene_expression',
	platform='Illumina HiSeq 2000', description='''
<p>Barplot of the expression of a single gene in multiple individual samples. 58 samples were sequenced using the whole transcriptome protocol of the Illumina HiSeq 2000 platform. Full 90 bp reads were aligned against the GRCh37 reference genome. Reads that aligned to exons annotated in Ensembl 68 were tallied for each gene. Gene expression values were normalized using median-of-ratios normalization.</p>
''')





# Prostate cancer - ExpO
dataset('Prostate cancer > IGC expO',
	title='Expression Project for Oncology',
	publication_date='23-12-2004',
	geo_series_id='GSE2109',
	platforms='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 Plus 2.0</a><br>')

expr_barplot('Prostate cancer > IGC expO > Gene expression barplot',
	category='Transcriptomic analysis', gene='ERG',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/intgen_expo',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 Plus 2.0</a><br>',
	description='Probe intensities were summarized into gene expressions using the RMA algorithm.')





# Friedlander et al. (2010)
dataset('Prostate cancer > Friedlander et al. (2012)',
	title='Friedlander et al. - Common structural and epigenetic changes in the genome of castration-resistant prostate cancer',
	publication_date='01-02-2012',
	pubmed_id='22158653',
	stats='Total samples: 15<br>Total patients: 14<br>',
	platforms='Agilent HG CGH 244A<br>Illumina Human Methylation 27<br>')

igv_track('Prostate cancer > Friedlander et al. (2012) > ' + 
	'All tumors, paired with pooled gDNA from healthy anonymous donors',
	category='Copy number tracks',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/friedlander12_cgh_tracks/tumor_vs_ref.seg.gz')

	





# Liu et al. (2009)
dataset('Prostate cancer > Liu et al. (2009)',
	title='Liu et al. - Copy number analysis indicates monoclonal origin of lethal metastatic prostate cancer',
	publication_date='12-04-2009',
	pubmed_id='19363497',
	stats='Total samples: 74<br>Total patients: 14<br>',
	platforms='Affymetrix Genome Wide SNP 6.0',
	content='''
Liu et al. show through a high-resolution genome-wide single nucleotide polymorphism and copy number survey that most, if not all, metastatic prostate cancers have monoclonal origins and maintain a unique signature copy number pattern of the parent cancer cell while also accumulating a variable number of separate subclonally sustained changes. They find no relationship between anatomic site of metastasis and genomic copy number change pattern. Taken together with past animal and cytogenetic studies of metastasis and recent single-locus genetic data in prostate and other metastatic cancers, these data indicate that despite common genomic heterogeneity in primary cancers, most metastatic cancers arise from a single precursor cancer cell.
''')

feature_matrix('Prostate cancer > Liu et al. (2009) > Feature matrix',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/bova_cna',
	platform='Affymetrix Genome Wide SNP 6.0', sample_labelsize='200')
	
igv_track('Prostate cancer > Liu et al. (2009) > ' + 
	'All tumors, paired with averaged non-cancerous samples',
	category='Copy number tracks',
	data='http://betastasis.cs.tut.fi/betadata/bova_cgh_tracks/bova_cna.seg')







# LuCaP xenografts
dataset('Prostate cancer > LuCaP xenografts',
	title='LuCaP xenografts', hidden=True,
	platforms='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 Plus 2.0</a><br>Agilent Human miRNA 8x15K v2<br>Agilent HG CGH 244A',
	content='''
<p>Gene expression for 22 LuCaP xenografts was measured using Affymetrix HG U133 Plus 2.0 microarrays. MicroRNA expression for 15 LuCaP xenografts was measured using Agilent Human miRNA 8x15K v2 microarrays.</p>
<p>Copy number profiles were measured for 19 LuCaP xenografts using Agilent HG CGH 244A microarrays. Xenografts were hybridized to the Cy5 channel and paired reference DNA to the Cy3 channel.</p>
''')

expr_barplot('Prostate cancer > LuCaP xenografts > Gene and miRNA expression barplot',
	category='Visualizations', gene='ERG',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/lucaps_expr',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 Plus 2.0</a><br>Agilent Human miRNA 8x15K v2',
	description='''
<p>Gene expression was measured for 22 LuCaP xenografts using Affymetrix HG U133 Plus 2.0 microarrays. MicroRNA expression was measured for 15 LuCaP xenografts using Agilent Human miRNA 8x15K v2 microarrays. Gene expression values were quantile normalized across all samples.</p>
''')

gene_boxplot('Prostate cancer > LuCaP xenografts > Gene and miRNA expression boxplot (fusion & AR)',
	category='Visualizations', gene='ERG',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/lucaps_expr',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/lucaps_groups/groups.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 Plus 2.0</a><br>Agilent Human miRNA 8x15K v2',
	description='''
<p>Gene expression was measured for 22 LuCaP xenografts using Affymetrix HG U133 Plus 2.0 microarrays. MicroRNA expression was measured for 15 LuCaP xenografts using Agilent Human miRNA 8x15K v2 microarrays. Gene expression values were quantile normalized across all samples.</p>
''')

feature_matrix('Prostate cancer > LuCaP xenografts > Feature matrix',
	category='Visualizations',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/lucaps_cna',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/lucaps_groups/groups.json',
	platform='Agilent HG CGH 244A', sample_labelsize='100',
	default_features='AR (CNA),TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),ERBB2 (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)')

igv_track('Prostate cancer > LuCaP xenografts > All xenografts, paired with Cy3 channel reference DNA',
	category='Copy number tracks',
	data='http://betastasis.cs.tut.fi/betadata/lucap_cgh_tracks/lucap_xenografts.seg')











# Massie et al. (2011)
dataset('Prostate cancer > Massie et al. (2011)',
	title='Massie et al. - The androgen receptor fuels prostate cancer by regulating central metabolism and biosynthesis',
	publication_date='20-05-2011',
	pubmed_id='21602788',
	geo_series_id='GSE18684',
	platforms='Illumina Human 6 v2 BeadArray',
	content='''
<p>Detailed analysis of androgen regulated gene expression in the LNCaP prostate cancer cell line. Since androgens and the AR are known to be important for prostate cancer cell proliferation and invasion we aimed to identify androgen receptor (AR) regulated genes by combining this detailed Illumina beadarray study of androgen regulated gene expression with AR ChIP-sequencing data.</p>

<p>LNCaP cells were grown in RPMI medium supplemented with 10% charcoal dextran stripped (steroid depleted) FBS for 72h, before treatment with 1nM R1881 or vehicle control. Total RNA was harvested every 30min for 4h and then every hour up to 24h using trizol (Sigma), quantified using a Nanodrop spectrophotometer (ND-1000). RNA samples were prepared for analysis on Illumina Human 6 v2 BeadArrays and Biotrove Realtime PCR panels according to the manufacturers protocols.</p>
''')

timeseries(
	'Prostate cancer > Massie et al. (2011) > Gene expression time series',
	category='Gene expression', gene='AR',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/androgen_exposure',
	control_group='LNCaP no treatment',
	platform='Illumina Human 6 v2 BeadArray',
	description='''
<p>Time series of the relative expression of a gene after exposure to 0.1% ethanol or 1nM R1881. Dark colored lines indicate averaged expression across replicates, relative to untreated LNCaP cells. The lightly shaded areas indicate the standard error of the mean.</p>

<p>From the paper: Forty-eight total RNA samples were harvested from LNCaP cells grown for 72 h in steroid depleted medium (RPMI supplemented with 10% charcoal dextran-stripped FBS). These comprised 3 time zero samples; 10 vehicle (ethanol) control samples taken at 2, 4, 8, 12 and 24 h in duplicate; 36 androgen (R1881)-treated samples taken every 30 min for 4 h then every hour until 24 h following treatment (with replicates at 1, 2, 4, 8, 12, 16, 20 and 24 h).</p>
''')











# Prostate cancer - Prensner et al. (2011)
dataset('Prostate cancer > Prensner et al. (2011)',
	title='Prensner et al. (2011) - Transcriptome sequencing across a prostate cancer cohort identifies PCAT-1, an unannotated lincRNA implicated in disease progression',
	publication_date='08-11-2010',
	geo_series_id='GSE25183',
	platforms='Illumina Genome Analyzer<br>Illumina Genome Analyzer II')

expr_barplot('Prostate cancer > Prensner et al. (2011) > Gene expression barplot',
	category='Transcriptomic analysis', gene='ERG',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/prensner11/gene_expression',
	platform='Illumina Genome Analyzer<br>Illumina Genome Analyzer II',
	description='Reads were aligned with Tophat against the GRCh37 genome assembly. Aligned read counts were then calculated against genes from Ensembl release 67. All expression values were quantile normalized.')










		
# Prostate cancer - Taylor et al. (2010)
dataset('Prostate cancer > Taylor et al. (2010)',
	title='Taylor et al. (2010) - Integrative Genomic Profiling of Human Prostate Cancer',
	publication_date='24-06-2010',
	pubmed_id='20579941',
	geo_series_id='GSE21032',
	stats='Total samples: 218<br>Total patients: 218<br>',
	platforms='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a><br>Agilent HG CGH 244A<br>Agilent miRNA 8x15K v2.0<br>',
	content='''
Taylor et al. concordantly assessed DNA copy number, exon expression and miRNA expression in 218 prostate cancer tumors. They identified the nuclear receptor coactivator NCOA2 as an oncogene in 11% of tumors. Additionally, the androgen-driven TMPRSS2-ERG fusion was associated with a previously unrecognized, prostate-specific deletion at chromosome 3p14 that implicates FOXP1, RYBP, and SHQ1 as potential cooperative tumor suppressors.
''')

gene_boxplot('Prostate cancer > Taylor et al. (2010) > Gene expression boxplot',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a><br>Agilent Human miRNA 8x15K v2')

kaplan_meier('Prostate cancer > Taylor et al. (2010) > Kaplan-Meier survival plot',
	category='Gene-level views',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_groups/grade.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a><br>Agilent Human miRNA 8x15K v2',
	recurrence='yes', feature='AR')

gene_scatter('Prostate cancer > Taylor et al. (2010) > Two-gene scatterplot',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression',
	gene_x='ETV1', gene_y='ERG',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression/groups/groups_new.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

feature_matrix('Prostate cancer > Taylor et al. (2010) > Feature matrix',
	category='Gene-level views',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_pca_features',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_groups/grade_trimmed.json',
	platform='Agilent HG CGH 244A<br>Sanger sequencing<br>iPlex Sequenom',
	sample_labelsize='70',
	default_features='AR (CNA),AR (MUT),TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)')

parallel_coordinates('Prostate cancer > Taylor et al. (2010) > Parallel coordinates',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_groups/grade.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a><br>Agilent Human miRNA 8x15K v2',
	features='["hsa-miR-143", "RB1", "E2F1", "AR", "MYC"]')

geneset_correlation('Prostate cancer > Taylor et al. (2010) > Gene set correlation',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression',
	scatter_url='/prostate_cancer/taylor_et_al/two-gene_scatterplot',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression/groups/groups_new.json',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>', ref_gene='AR',
	target_genes='AURKA,KLK3,ERG,ETV1,ETV4,PTEN,TP53,NCOA2')

alt_splicing(
	'Prostate cancer > Taylor et al. (2010) > Alternative splicing boxplot (gene)',
	category='Alternative splicing', gene='AR',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_gene_expression/groups/groups_new.json',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/taylor_splice_new',
	default_test_group='Primary (untreated)',
	default_ref_group='Normal (BPH)',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

alt_splicing_exons('Prostate cancer > Taylor et al. (2010) > Alternative splicing boxplot (exon)',
	category='Alternative splicing',
	gene='AR', data='http://betastasis.cs.tut.fi/betadata/taylor_splice_new',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

alt_splicing_survival(
	'Prostate cancer > Taylor et al. (2010) > Alternative splicing survival association',
	category='Alternative splicing', gene='AR', recurrence='yes',
	data='http://betastasis.cs.tut.fi/betadata/taylor_splice_new',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_human_exon_10_st/">Affymetrix Human Exon 1.0 ST</a>')

igv_track('Prostate cancer > Taylor et al. (2010) > All tumors, paired with adjacent normal or universal reference DNA',
	category='Copy number tracks',
	data='http://betastasis.cs.tut.fi/betadata/taylor_pca_all.seg')








# Prostate cancer - Tamura (2007)
dataset('Prostate cancer > Tamura (2007)',
	title='Tamura et al. (2007) - Molecular features of hormone-refractory prostate cancer cells by genome-wide gene-expression profiles',
	publication_date='01-06-2007',
	geo_series_id='GSE6811',
	platforms='Custom cDNA microarray<br>')

expr_barplot('Prostate cancer > Tamura (2007) > Gene expression barplot',
	category='Transcriptomic analysis', gene='ETV4',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tamura07',
	platform='Custom cDNA microarray<br>',
	description='Tumor samples were hybridized onto a custom cDNA microarray. Expression values come directly from the GEO series matrix.')











# Prostate cancer - TCGA PRAD
dataset('Prostate cancer > TCGA PRAD',
	title='TCGA - Prostatic adenocarcinoma (PRAD)',
	stats='Total samples: 130<br>Tumor samples: 65<br>Normal blood samples: 55<br>Adjacent normal samples: 10<br>Total patients: 65<br>',
	platforms='Affymetrix Genome Wide SNP 6.0',
	content='''
The TCGA prostate cancer project has currently produced measurements for 130 samples using the Affymetrix GW SNP 6.0 array.
''')

expr_barplot('Prostate cancer > TCGA PRAD > Gene expression barplot (RNA-seq)',
	category='Transcriptome sequencing', gene='AR',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tcga_prad/rnaseq_gene_expression',
	platform='Illumina Genome Analyzer II', description='''
<p>Quantile normalized expression values for genes annotated in Ensembl release 67.</p>
''')

feature_matrix('Prostate cancer > TCGA PRAD > Feature matrix',
	category='Gene-level views', data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tcga_prad/features',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tcga_prad/groups/groups.json', 
	platform='Affymetrix Genome Wide SNP 6.0', sample_labelsize='100',
	default_features='AR (CNA),TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)')
	
igv_track('Prostate cancer > TCGA PRAD > ' +
	'All tumors, paired with adjacent normal or universal reference',
	category='Copy number tracks',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tcga_prad/tcga_pca_snp6.seg')






# Prostate cancer - Tampere PC
dataset('Prostate cancer > Tampere PC',
	title='Tampere PC sequencing', hidden=True,
	stats='RNA-seq samples: 53<br>sRNA-seq samples: 45<br>DNA-seq samples: 40<br>MeDIP-seq samples: 51<br><br><b>Links</b><br><a href="http://sumatranlehvi.cs.tut.fi/RE/">Regulome Explorer</a><br>',
	platforms='Illumina HiSeq 2000<br>Solexa (sRNA)<br>',
	content='''
<p>This dataset contains multimodal high throughput sequencing data from 28 untreated prostate cancers (PC) and 13 castration resistant prostate cancers (CRPC). 12 benign prostatic hyperplasias (BPH) were also sequenced as controls. The samples were sequenced using transcriptome (RNA-seq), small RNA (sRNA-seq), low-coverage whole genome (DNA-seq), and methylated DNA sequencing (MeDIP-seq). A <b><a href="http://sumatranlehvi.cs.tut.fi/RE/">Regulome Explorer visualization</a></b> of the dataset is now available.</p>
''')

expr_barplot('Prostate cancer > Tampere PC > Gene expression barplot',
	category='Whole transcriptome analysis', gene='AR',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/rnaseq_gene_expression',
	platform='Illumina HiSeq 2000', description='''
<p>53 samples were sequenced using the whole transcriptome protocol of the Illumina HiSeq 2000 platform. Full 90 bp reads were aligned against NCBI RefSeq 38 transcripts sequences. 2 nucleotide mismatches were allowed in alignments. 70-80% of all reads aligned against annotated transcript sequences. Gene expression values were quantile normalized across all 53 samples.
</p>
''')

expr_barplot('Prostate cancer > Tampere PC > Gene expression barplot (Ensembl 68)',
	category='Whole transcriptome analysis', gene='AR',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/gene_expression_ensembl',
	platform='Illumina HiSeq 2000', description='''
<p>53 samples were sequenced using the whole transcriptome protocol of the Illumina HiSeq 2000 platform. Full 90 bp reads were aligned against the GRCh37 genome with Tophat. Gene expression values were calculated using bedtools based on composite gene exon models. Gene expression values were normalized using median-of-ratios normalization.</p>
''')

gene_boxplot('Prostate cancer > Tampere PC > Gene expression boxplot',
	category='Whole transcriptome analysis', gene='AR',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/rnaseq_gene_expression',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/groups/tumor_grade_groups.json',
	platform='Illumina HiSeq 2000',
	description='''
<p>53 prostate cancer samples were sequenced using the whole transcriptome protocol of the Illumina HiSeq 2000 platform. Full 90 bp reads were aligned against NCBI RefSeq 38 transcripts sequences. 2 nucleotide mismatches were allowed in alignments. 70-80% of all reads aligned against annotated transcript sequences.
''')

#exon_expr_barplot('Prostate cancer > Tampere PC > Exon expression barplot',
#	category='Whole transcriptome analysis', gene='AR', exon='1',
#	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/visakorpi_bgi_exon_expr',
#	platform='Illumina HiSeq 2000', description='''
#<p>53 samples were sequenced using the whole transcriptome protocol of the Illumina HiSeq 2000 platform. Full 90 bp reads were aligned against NCBI RefSeq 38 transcripts sequences. 2 nucleotide mismatches were allowed in the alignments. 70-80% of all reads aligned against annotated transcript sequences. The expression of an exon was calculated by summing all reads that overlapped with the exon's sequence within a transcript.Total read counts were quantile normalized across samples.</p>
#''')

gene_scatter('Prostate cancer > Tampere PC > Two-gene scatterplot',
	category='Whole transcriptome analysis', gene_x='ERG', gene_y='ETV1',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/rnaseq_gene_expression',
	groups='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/groups/tumor_grade_groups.json',
	platform='Illumina HiSeq 2000')

#graph('Prostate cancer > Tampere PC > Outlier association graph',
#	category='Whole transcriptome analysis',
#	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca_associations.json')

expr_barplot('Prostate cancer > Tampere PC > MicroRNA expression barplot',
	category='Small RNA analysis', gene='hsa-miR-143-5p',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/mirna_expression',
	feature_type='miRNA',
	platform='Solexa')

#expr_barplot('Prostate cancer > Tampere PC > MicroRNA expression barplot (ProspeR)',
#	category='Small RNA analysis', gene='hsa-miR-21',
#	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/visakorpi_mirna_expr',
#	platform='Solexa')

igv_track('Prostate cancer > Tampere PC > Copy number tracks (IGV)',
	category='Whole genome sequencing',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/tampere_pca_cnv_seq_cbs.seg')

circvis('Prostate cancer > Tampere PC > Chromosomal rearrangements',
	category='Whole genome sequencing', gene='', default='PC_6864',
	title='Chromosomal rearrangements',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/tampere_pca/structural_variants',
	platform='Illumina HiSeq 2000')

igv_track('Prostate cancer > Tampere PC > Methylation tracks (IGV)',
	category='Methylated DNA immunoprecipitation sequencing',
	data='http://viherjora.cs.tut.fi/PCa/coverage_tracks/meth_coverages_corrected.tdf')


#aberrated_scatter('Prostate cancer > Tampere PC > Scatterplot of most aberrated genes (HNPCa vs BPH)',
#	category='Integrative analysis',
#	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/top_genes_hnpca_vs_bph.json',
#	platform='Illumina HiSeq 2000')

aberrated_scatter('Prostate cancer > Tampere PC > Scatterplot of most aberrated genes (CRPC vs HNPCa)',
	category='Integrative analysis',
	data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/top_genes_crpc_vs_pca.json',
	platform='Illumina HiSeq 2000')


#mutation_list('Prostate cancer > Tampere PC > Mutation list',
#	category='Mutation analysis', data='http://betastasis.cs.tut.fi/betadata/prostate_cancer/')





###########
# TISSUES #
###########

# Affymetrix tissue dataset
dataset('Tissues > Affymetrix tissue dataset',
	title='Affymetrix tissue dataset')

expr_barplot('Tissues > Affymetrix tissue dataset > ' + 
	'Gene expression barplot (Affymetrix HG U133 Plus 2.0',
	category='Gene expression', gene='TP53',
	data='http://betastasis.cs.tut.fi/betadata/tissues/affy_u133',
	platform='<a href="/platforms/cdna_microarrays/affymetrix_hg_u133_plus_20/">Affymetrix HG U133 Plus 2.0</a>')

# Applied Biosystems tissue microRNA expression dataset
dataset('Tissues > ABI tissue microRNA expression',
	title='ABI tissue microRNA expression')

expr_barplot('Tissues > ABI tissue microRNA expression > ' + 
	'MicroRNA expression barplot',
	category='MicroRNA expression', gene='hsa-miR-21',
	data='http://betastasis.cs.tut.fi/betadata/tissues/abi_mirna_tissue_expr',
	platform='TaqMan microRNA assays')
	






dataset('Tissues > Cancer Cell Line Encyclopedia',
	title='Cancer Cell Line Encyclopedia',
	stats='Total samples: 967<br>',
	platforms='Affymetrix HG U133 Plus 2.0<br>',
	content='''
<p>Expression profiling of cancer cell lines.</p>
''')

expr_barplot('Tissues > Cancer Cell Line Encyclopedia > Gene expression barplot',
	category='Transcriptomic analysis', gene='TP53',
	data='http://betastasis.cs.tut.fi/betadata/tissues/ccle/gene_expression',
	groups='http://betastasis.cs.tut.fi/betadata/tissues/ccle/groups/tissue_groups.json',
	platform='Affymetrix HG U133 Plus 2.0', description='''
<p>Gene expression data was downloaded directly from the Cancer Cell Line Encyclopedia website. According to the website, raw Affymetrix CEL files were converted to a single value for each probeset using RMA and quantile normalization. The probeset used was ENTREZG v15 from the BrainArray project.</p>''')









#############
# PLATFORMS #
#############

# cDNA microarrays
dataset('Platforms > cDNA microarrays',
	title='cDNA microarrays', stats='Total arrays: 7<br>')

cdna_uarray_probes('Platforms > cDNA microarrays > Affymetrix HG U133 Plus 2.0',
	data='http://betastasis.cs.tut.fi/betadata/platforms/affymetrix_hg_u133_plus_20')
	
cdna_uarray_probes('Platforms > cDNA microarrays > Affymetrix HG U133A',
	data='http://betastasis.cs.tut.fi/betadata/platforms/affymetrix_hg_u133a')

cdna_uarray_probes('Platforms > cDNA microarrays > Affymetrix HT HG U133A',
	data='http://betastasis.cs.tut.fi/betadata/platforms/affymetrix_ht_hg_u133a')

cdna_uarray_probes(
	'Platforms > cDNA microarrays > Affymetrix Human Exon 1.0 ST',
	data='http://betastasis.cs.tut.fi/betadata/platforms/affymetrix_human_exon_10_st')

cdna_uarray_probes('Platforms > cDNA microarrays > Agilent 244K TCGA custom #1',
	data='http://betastasis.cs.tut.fi/betadata/platforms/agilent_244k_tcga_custom_1')

cdna_uarray_probes('Platforms > cDNA microarrays > Agilent 244K TCGA custom #2',
	data='http://betastasis.cs.tut.fi/betadata/platforms/agilent_244k_tcga_custom_2')

cdna_uarray_probes('Platforms > cDNA microarrays > Agilent 244K TCGA custom #3',
	data='http://betastasis.cs.tut.fi/betadata/platforms/agilent_244k_tcga_custom_3')








#########
# TOOLS #
#########

circvis('Tools > Circvis',
	category='Circvis', default='', data='')





##########################################
# INITIALIZATION OF THE WSGI APPLICATION #
##########################################

bottle.debug(True)
bottle.TEMPLATE_PATH = [ 'templates' ]

def wsgi_app():
	return bottle.default_app()
	
if __name__ == '__main__':
	bottle.run(app=wsgi_app(), host='localhost', port=8080)

