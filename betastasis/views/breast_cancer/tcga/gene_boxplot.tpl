%scripts = [ 'ui', 'gene_boxplot' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>TCGA - Breast cancer</h2>

%include gene_boxplot gene='TP53', data_root='/betadata/breast_cancer/tcga_expr', groups='/betadata/breast_cancer/tcga_groups/t_stage_boxplot.json', platform='<a href="/microarrays/agilent_244k_tcga_custom_3">Agilent 244K TCGA custom 3</a>'


