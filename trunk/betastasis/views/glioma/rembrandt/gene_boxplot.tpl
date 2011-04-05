%scripts = [ 'ui', 'gene_boxplot' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>Rembrandt</h2>

%include gene_boxplot gene='TP53', data_root='/betadata/rembrandt_foo', platform='<a href="/microarrays/affy_hg_u133_plus_2">Affymetrix HG U133 v2.0 Plus</a>'

