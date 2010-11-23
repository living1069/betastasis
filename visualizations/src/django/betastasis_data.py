from betastasis_tools import attrgetter

from django.utils.datastructures import SortedDict

DATA_SET_NAME = 'name'
DESCRIPTION = 'description'
CANCER_TYPE = 'cancer_type'
DATA_SET_ROOT = 'data_set_root'
DEFAULT_GENE = 'default_gene'
EXPERIMENT_PLATFORMS = 'experiment_platforms'
MENU_TEXT = 'menu_text'
FEATURE_LIST = 'feature_list'

PLATFORM_NAME = 'platform_name'
EXPRESSION_PATH = 'expression_path'
SURVIVAL_PATH = 'survival_path'

# Cancer types
# ============
cancer_types = {'glioma': 'Glioma',
                'prostate': 'Human Prostate Cancer',
                'ovarian': 'Ovarian Cancer'}

# Experiment platforms
# ====================
class Probeset:
    def __init__(self, platform_id, transcriptome_id, path):
        self.platform_id = platform_id
        self.transcriptome_id = transcriptome_id
        self.path = path

transcriptomes = {'ncbi_refseq_38': 'NCBI RefSeq 38'}

probesets = [Probeset('affy_ht_hg_u133a', 'ncbi_refseq_38', 'affy_ht_hg_u133a'),
             Probeset('affy_huex_10_st', 'ncbi_refseq_38', 'affy_huex'),
             Probeset('agilent_244k_tcga_custom_1', 'ncbi_refseq_38', 'agilent_244k_tcga_custom_1'),
             Probeset('agilent_244k_tcga_custom_2', 'ncbi_refseq_38', 'agilent_244k_tcga_custom_2'),
             Probeset('agilent_244k_tcga_custom_3', 'ncbi_refseq_38', 'agilent_244k_tcga_custom_3')]

platforms = {'affy_ht_hg_u133a': 'Affymetrix HT HG U133A',
             'affy_huex_10_st': 'Affymetrix Human Exon 1.0 ST',
             'affy_hg_u133_v20_plus': 'Affymetrix HG U133 v2.0 Plus',
             'agilent_244k_tcga_custom_1': 'Agilent 244K TCGA custom 1',
             'agilent_244k_tcga_custom_2': 'Agilent 244K TCGA custom 2',
             'agilent_244k_tcga_custom_3': 'Agilent 244K TCGA custom 3',
             'agilent_hg_cgh_244a': 'Agilent HG CGH 244A'}

# Data sources
# ============
data_sources = {'tcga_gbm': {DATA_SET_NAME: 'The Cancer Genome Atlas - Glioblastoma Multiforme',
                             DESCRIPTION: 'The TCGA GBM dataset consists of primary tumor samples from roughly 400 patients. The dataset contains data from gene expression, CGH and methylation microarrays.',
                             MENU_TEXT: 'TCGA GBM'},
                'rembrandt': {DATA_SET_NAME: 'REMBRANDT',
                              DESCRIPTION: 'REpository for Molecular BRAin Neoplasia DaTa (REMBRANDT) is a robust bioinformatics knowledgebase framework that leverages data warehousing technology to host and integrate clinical and functional genomics data from clinical trials involving patients suffering from gliomas.',
                              MENU_TEXT: 'REMBRANDT'},
                'taylor': {DATA_SET_NAME: 'Taylor et al. - Cancer Cell Volume 18, Issue 1, 13 July 2010, Pages 11-22',
                           DESCRIPTION: 'Integrative Genomic Profiling of Human Prostate Cancer',
                           MENU_TEXT: 'Taylor et al.'}}

# Paths to data set feature listings
# ==================================
data_set_features = { 'tcga_gbm': 'gbm_cna_features.json',
                      'taylor': 'prostate_cna_features.json' }

class DataSet:
    def __init__(self, cancer_type, data_type, data_source, platform, path, transcriptome):
        self.cancer_type = cancer_type
        self.data_type = data_type
        self.data_source = data_source
        self.platform = platform
        self.path = path
        self.transcriptome = transcriptome

data_sets = [DataSet('glioma', 'expression', 'tcga_gbm', 'affy_ht_hg_u133a', 'tcga_gbm_u133a_expr', 'ncbi_refseq_38'),
             DataSet('glioma', 'feature_matrix', 'tcga_gbm', 'agilent_hg_cgh_244a', 'gbm_cna_features.json', 'ncbi_refseq_38'),
             DataSet('glioma', 'survival', 'tcga_gbm', 'affy_ht_hg_u133a', 'tcga_gbm_u133_survival', 'ncbi_refseq_38'),
             DataSet('glioma', 'survival', 'tcga_gbm', 'affy_huex_10_st', 'tcga_gbm_huex_survival', 'ncbi_refseq_38'),
             DataSet('glioma', 'expression', 'rembrandt', 'affy_hg_u133_v20_plus', 'rembrandt_expr', 'ncbi_refseq_38'),
             DataSet('glioma', 'survival', 'rembrandt', 'affy_hg_u133_v20_plus', 'rembrandt_survival', 'ncbi_refseq_38'),
             DataSet('prostate', 'expression', 'taylor', 'affy_huex_10_st', 'taylor_prostate_expr', 'ncbi_refseq_38'),
             DataSet('prostate', 'feature_matrix', 'taylor', 'agilent_hg_cgh_244a', 'prostate_cna_features.json', 'ncbi_refseq_38')]

def fill_data_source_meta(p_data_sets):
    filled = []
    for ds in p_data_sets:
        row = {}
        row['cancer_type'] = {'id': ds.cancer_type, DESCRIPTION: cancer_types[ds.cancer_type] }
        row['data_type'] = ds.data_type
        row['data_source'] = {'id': ds.data_source, 'meta': data_sources[ds.data_source]}
        row['platform'] = {'id': ds.platform, 'name': platforms[ds.platform]}
        row['path'] = ds.path
        row['transcriptome'] = {'id': ds.transcriptome, 'name': transcriptomes[ds.transcriptome]}
        
        # Does probeset information exist for this platform and transcriptome?
        if get_probeset_path_if_exists(ds.platform, ds.transcriptome) is not None:
            row['platform']['probeset_exists'] = True
        
        filled.append(row)

    return filled

def get_data_sets_by_cancer_type():
    ds = sorted(data_sets, key=attrgetter('cancer_type', 'data_source', 'data_type'))
    return fill_data_source_meta(ds)

def make_current_cancer_dict(ct):
    d = {'id': ct, 'name': cancer_types[ct]}
    return d

def get_probeset_path_if_exists(platform_id, transcriptome_id):
    for ps in probesets:
        if (ps.platform_id == platform_id and
            ps.transcriptome_id == transcriptome_id):
            return ps.path
    
    return None

def get_probesets_by_platform():
    temp = sorted(probesets, key=attrgetter('platform_id'))
    filled = []
    for ps in temp:
        row = {}
        row['platform'] = {'id': ps.platform_id, 'name': platforms[ps.platform_id]}
        row['transcriptome'] = {'id': ps.transcriptome_id, 'name': transcriptomes[ps.transcriptome_id]}
        filled.append(row)
    
    return filled

def get_data_path_if_exists(data_source_id, platform_id, data_type):
    for ds in data_sets:
        if (ds.data_source == data_source_id and
            ds.platform == platform_id and
            ds.data_type == data_type):
            return ds.path
    
    return None
