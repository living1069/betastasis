
<h3>Gene expression bar plot</h3>
<div id="fig"></div>
</div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
Gene: <input type="text" id="GeneSelect" value="" /><br/><br/>

<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br>
</p>
</div>

</div>
	  
<script type="text/javascript">
var data_root = '{{data_root}}';

$(document).ready(function() {
	
	var data = {gene: ''};
	var default_gene = gup('gene');
	if (default_gene.length == 0)
		default_gene = '{{gene if 'gene' in locals() else 'TP53'}}';
	$('#GeneSelect').val(default_gene);
	
	var data_downloaded = function() {
		if (!('sample_ids' in data && 'expr' in data))
			return;
			
		vis.gene_name = data.gene;
		vis.set_data(data.expr, data.sample_ids);
	}

	var load_gene_data = function(gene_name) {
		if (gene_name == data.gene) return;
		$.getJSON(data_root + '/expr/' + gene_name[0].toLowerCase() + '/'
			+ gene_name + '.json', function(d) {
			data.expr = d['data'];
			data.gene = gene_name;
			data_downloaded();
		});
	}

	var vis = new Barplot('fig');

    load_gene_data(default_gene);
	
    $.getJSON(data_root + '/expr/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			load_gene_data(val);
		});
	});
	
	$.getJSON(data_root + '/clinical/sample_id.json', function(d) {
		data.sample_ids = d['data']; data_downloaded();
	});

});
</script>

