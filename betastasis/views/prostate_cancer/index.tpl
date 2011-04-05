%rebase base nav_tree=nav_tree	

<link rel="stylesheet" type="text/css" href="/style/tooltip.css" />
<script type="text/javascript" language="javascript" src="/js/pathway_atlas.js"></script>

<div id="main">

<h2>Prostate cancer</h2>
		
<p>Prostate cancer is a form of cancer that develops in the prostate, a gland in the male reproductive system. Most prostate cancers are slow growing; however, there are cases of aggressive prostate cancers. The cancer cells may metastasize from the prostate to other parts of the body, particularly the bones and lymph nodes. Prostate cancer usually progresses from adenocarcinoma in situ to adenocarcinoma and then metastatic disease.</p>

<h3>Available datasets</h3><ul>
	<li><a href="/prostate_cancer/bova/">Bova et al.</a></li>
	<li><a href="/prostate_cancer/lucap_xenografts/">LuCaP xenografts</a></li>
	<li><a href="/prostate_cancer/taylor/">Taylor et al.</a></li>
</ul>

<h3>Pathways</h3>
<object id="test_svg" onload="$.getJSON('/betadata/cancer_pathways/pca_annotations.json', function(d) { annotations_loaded(d); });" data="/betadata/cancer_pathways/prostate_cancer.svgz" type="image/svg+xml" width="900" height="1500"></object>

</div>

%include news

