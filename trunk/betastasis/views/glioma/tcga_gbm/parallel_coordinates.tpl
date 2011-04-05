%scripts = [ 'ui', 'parallel_coordinates', 'feature_selector' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>TCGA - Glioblastoma multiforme</h2>

%include parallel_coordinates data_root='/betadata/tcga_gbm_u133', platform='<a href="/microarrays/affy_ht_hg_u133a">Affymetrix HT HG U133A</a>', features='["EGFR", "PTEN", "TP53", "PDGFRA", "IDH1"]'

