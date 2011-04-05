%scripts = [ 'ui', 'parallel_coordinates', 'feature_selector' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>Taylor et al. - Integrative Genomic Profiling of Human Prostate Cancer</h2>

%include parallel_coordinates data_root='/betadata/taylor_foo', groups='/betadata/prostate_cancer/taylor_groups/grade.json', platform='<a href="/microarrays/affy_huex">Affymetrix Human Exon 1.0 ST</a><br>Agilent Human miRNA 8x15K v2', features='["hsa-miR-143", "RB1", "E2F1", "AR", "MYC"]'

