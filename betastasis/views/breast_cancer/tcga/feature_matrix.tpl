%scripts = [ 'ui', 'feature_matrix' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>TCGA - Breast cancer</h2>

%include feature_matrix data_root='/betadata/breast_cancer/tcga_cna', groups='/betadata/breast_cancer/tcga_groups/t_stage.json', platform='Affymetrix GW SNP 6.0', sample_labelsize='70', default_features='TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),ERBB2 (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)'

