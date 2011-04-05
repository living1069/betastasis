<h3>Probeset design</h3>
<div id="fig"></div>
</div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Platform: {{platform if 'platform' in locals() else '-'}}<br>
Transcriptome: NCBI RefSeq 38<br><br>
Gene: <input type="text" id="GeneSelect" value="TP53" />
</p>
</div>
</div>

<script type="text/javascript">
var genelist;
var data_root = '{{data_root}}';
	
$(document).ready(function() {
	var vis = new ProbesetVis('fig');
	
	var genelist_downloaded = function(d) {
		genelist = d['genes'];
		$('#GeneSelect').omnicomplete(genelist, function(val) {
			$.getJSON(data_root + '/' + val[0].toLowerCase() + '/'
				+ val + '.json', function(d) { vis.setGene(d); });
		});
	}
	
	$.getJSON(data_root + '/t/TP53.json',
		function(d) { vis.setGene(d); });
	$.getJSON(data_root + '/genelist.json', genelist_downloaded);
});
</script>

