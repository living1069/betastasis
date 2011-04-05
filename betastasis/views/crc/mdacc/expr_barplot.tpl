%scripts = [ 'ui', 'barplot' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>MDACC - Colon cancer RNA-seq</h2>

%include gene_expr_barplot gene='TP53', data_root='/betadata/colorectal_cancer/mda_rnaseq_expr', platform='SOLiD 3'


