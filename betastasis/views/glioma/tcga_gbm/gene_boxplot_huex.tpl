%scripts = [ 'ui', 'gene_boxplot' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>TCGA - Glioblastoma multiforme</h2>

%include gene_boxplot gene='TP53', data_root='/betadata/tcga_gbm_huex', platform='<a href="/microarrays/affy_huex">Affy Human Exon 1.0 ST</a>'


