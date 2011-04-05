%scripts = [ 'ui', 'feature_matrix' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>Taylor et al. - Somatic mutations of the Parkinson's disease-associated gene PARK2 in glioblastoma and other human malignancies</h2>

%include feature_matrix data_root='/betadata/colorectal_cancer/taylor_cna', groups='/betadata/colorectal_cancer/taylor_groups/groups.json', platform='Agilent HG CGH 244A', sample_labelsize='70', default_features='TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)'

