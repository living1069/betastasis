%rebase base nav_tree=nav_tree

<link rel="stylesheet" type="text/css" href="/style/tooltip.css" />
<script type="text/javascript" language="javascript" src="/js/pathway_atlas.js"></script>

<div id="main">

<h2>Glioma</h2>

<p>A glioma is a type of tumor that arises from glial cells. Gliomas usually occur in the brain, but sometimes also in the spine.</p>
		
<h3>Available datasets</h3><ul>
	<li><a href="/glioma/rembrandt">Rembrandt</a></li>
	<li><a href="/glioma/tcga_gbm">TCGA GBM</a></li>
</ul>

<h3>Classification</h3>
<p>Gliomas are classified according to the type of cell from which the tumor originally arose:<ul>
<li>Astrocytoma: astrocytes</li>
<li>Oligodendroglioma: oligodendrocytes</li>
<li>Oligoastrocytoma: mixed appearance, features from both oligodendrocytes and astrocytes</li>
<li>Ependymomas: ependymal cells</li>
</ul></p>
		
<p>Gliomas are also classified by tumor grade into low grade (well-differentiated) and high grade (anaplastic) tumors.</p>

<h3>Pathways</h3>
<object id="test_svg" onload="$.getJSON('/betadata/cancer_pathways/gbm_annotations.json', function(d) { annotations_loaded(d); });" data="/betadata/cancer_pathways/gbm.svgz" type="image/svg+xml" width="900" height="1500"></object>

</div>

%include news

