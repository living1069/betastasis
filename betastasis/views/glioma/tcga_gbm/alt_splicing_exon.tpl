%scripts = [ 'ui', 'gene_boxplot' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>TCGA - Glioblastoma multiforme</h2>

%include exon_boxplot gene='TP53', data_root='/betadata/glioma/tcga_gbm_splice', groups='/betadata/glioma/tcga_gbm_groups/groups.json', platform='<a href="/microarrays/affy_huex">Affymetrix Human Exon 1.0 ST</a>'

</div>
