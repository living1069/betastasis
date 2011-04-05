%scripts = [ 'ui', 'feature_matrix' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<div id="main">
<h2>LuCaP xenografts</h2>

%include feature_matrix data_root='/betadata/lucaps_cna', platform='Agilent HG CGH 244A', sample_labelsize='100'

