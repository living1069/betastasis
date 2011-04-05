import bottle
import os
import os.path
from bottle import route, run, template, static_file, abort, redirect

bottle.debug(True)
bottle.TEMPLATE_PATH = [ 'templates', 'views' ]

hierarchy = {
	'breast_cancer': {
		'index': 'Breast cancer',
		'tcga': {
			'index': 'TCGA BRCA',
			'gene_boxplot': 'Gene boxplot',
			'gene_survival': 'Gene survival association',
			'parallel_coordinates': 'Parallel coordinates',
			'feature_matrix': 'Feature matrix',
		},
	},
	'crc': {
		'index': 'Colorectal cancer',
		'taylor': {
			'index': 'Taylor et al.',
			'feature_matrix': 'Feature matrix',
		},
		'mdacc': {
			'index': 'MDACC',
			'expr_barplot': 'Gene expression bar plot',
		},
	},
	'glioma': {
		'index': 'Glioma',
		'tcga_gbm': {
			'index': 'TCGA GBM',
			'feature_matrix': 'Feature matrix',
			'gene_boxplot_u133': 'Gene boxplot (Affymetrix HT HG U133A)',
			'gene_boxplot_huex': 'Gene boxplot (Affymetrix Human Exon 1.0 ST)',
			'parallel_coordinates': 'Parallel coordinates',
			'gene_survival_u133': 'Gene survival association (Affymetrix HT HG U133A)',
			'gene_survival_huex': 'Gene survival association (Affymetrix Human Exon 1.0 ST)',
			'alt_splicing_survival': 'Alt splicing survival association',
			'alt_splicing': 'Alt splicing boxplot (gene)',
			'alt_splicing_exon': 'Alt splicing boxplot (exon)',
			'alt_splicing_toplist': 'Alt splicing toplist',
		},
		'rembrandt': {
			'index': 'Rembrandt',
			'gene_boxplot': 'Gene boxplot',
			'gene_survival': 'Gene survival association',
			'parallel_coordinates': 'Parallel coordinates',
		},
	},
	'prostate_cancer': {
		'index': 'Prostate cancer',
		'taylor': {
			'index': 'Taylor et al',
			'feature_matrix': 'Feature matrix',
			'gene_boxplot': 'Gene boxplot',
			'parallel_coordinates': 'Parallel coordinates',
			'gene_survival': 'Gene survival association',
			'alt_splicing_survival': 'Alt splicing survival association',
			'alt_splicing': 'Alt splicing boxplot (gene)',
			'alt_splicing_exon': 'Alt splicing boxplot (exon)',
			'alt_splicing_toplist': 'Alt splicing toplist',
		},
		'bova': {
			'index': 'Bova et al',
			'feature_matrix': 'Feature matrix',
		},
		'lucap_xenografts': {
			'index': 'LuCaP xenografts',
			'feature_matrix': 'Feature matrix',
		},
	},
	'microarrays': {
		'index': 'Microarrays',
		'affy_hg_u133a': {'index': 'Affymetrix HG U133A' },
		'affy_hg_u133_plus_2': { 'index': 'Affymetrix HG U133 Plus 2.0' },
		'affy_ht_hg_u133a': { 'index': 'Affymetrix HT HG U133A' },
		'affy_huex': { 'index': 'Affymetrix Human Exon 1.0 ST' },
		'agilent_244k_tcga_custom_1': { 'index': 'Agilent 244K TCGA custom #1'},
		'agilent_244k_tcga_custom_2': { 'index': 'Agilent 244K TCGA custom #2'},
		'agilent_244k_tcga_custom_3': { 'index': 'Agilent 244K TCGA custom #3'},
	},
}

@route('/js/:filename#.+#')
def get_script(filename):
	return static_file(filename, root='js')

@route('/style/:filename#.+#')
def get_style(filename):
	return static_file(filename, root='style')

@route('/images/:filename#.+#')
def get_image(filename):
	return static_file(filename, root='images')

@route('/betadata/:filename#.+#')
def redirect_data(filename):
	return static_file(filename, root='/home/csbgroup/public_html/betadata')


@route('/')
def index():
	return template('index', nav_tree={'levels': []})

@route('/:page#.+#')
def get_page(page):
	levels = page.split('/')
	levels = [level for level in levels if len(level) >= 1]
	node = hierarchy
	
	for level in levels:
		if not level in node:
			abort(404, 'Sorry, that page does not exist.')
		node = node[level]
	
	nav_tree = hierarchy.copy()
	nav_tree['levels'] = levels
	
	if type(node) is dict:
		return template('/'.join(levels) + '/index', nav_tree=nav_tree)
	else:
		return template('/'.join(levels), nav_tree=nav_tree)
	
	#for level in levels[0:-1]:
	#	if not os.path.isdir('views' + curr_dir + '/' + level):
	#		abort(404, 'Sorry, that page does not exist.')
	#	curr_dir += '/' + level
	#	for fname in os.listdir('views/' + curr_dir)


run(host='localhost', port=8080)

