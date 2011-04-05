%scripts = [ 'ui', 'feature_matrix' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>Taylor et al. - Integrative Genomic Profiling of Human Prostate Cancer</h2>

%include feature_matrix data_root='/betadata/prostate_cancer/taylor_features', groups='/betadata/prostate_cancer/taylor_groups/grade_trimmed.json', platform='Agilent HG CGH 244A<br>Sanger sequencing<br>iPlex Sequenom', sample_labelsize='70', default_features='AR (CNA),AR (MUT),TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)'

