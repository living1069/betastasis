%scripts = [ 'ui', 'feature_matrix' ]
%rebase base nav_tree=nav_tree, scripts=scripts
	
<div id="main">
<h2>TCGA - Glioblastoma multiforme</h2>

%include feature_matrix data_root='/betadata/glioma/tcga_gbm_cna', groups='/betadata/glioma/tcga_gbm_groups/groups.json', platform='Agilent HG CGH 244A', default_features='EGFR (CNA),PDGFRA (CNA),CDK4 (CNA),CDK6 (CNA),CDKN2A (CNA),MDM2 (CNA),MDM4 (CNA),PIK3C2B (CNA),AKT1 (CNA),TP53 (CNA),PTEN (CNA),MYC (CNA),MET (CNA)'


