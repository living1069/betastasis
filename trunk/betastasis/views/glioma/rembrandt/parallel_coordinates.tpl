%scripts = [ 'ui', 'parallel_coordinates', 'feature_selector' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>Rembrandt</h2>

%include parallel_coordinates data_root='/betadata/rembrandt_foo', platform='<a href="/microarrays/affy_hg_u133_plus_2">Affymetrix HG U133 v2.0 Plus</a>', features='["EGFR", "PTEN", "TP53", "PDGFRA", "IDH1"]'

