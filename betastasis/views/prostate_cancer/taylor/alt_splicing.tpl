%scripts = [ 'ui', 'gene_boxplot' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>Taylor et al. - Integrative Genomic Profiling of Human Prostate Cancer</h2>

%include splice_boxplot gene='AR', data_root='/betadata/taylor_splice', default_test_group='Primary (untreated)', default_ref_group='Normal (BPH)',  platform='<a href="/microarrays/affy_huex">Affymetrix Human Exon 1.0 ST</a>'

</div>
