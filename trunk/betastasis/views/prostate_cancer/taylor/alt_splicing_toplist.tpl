%scripts = [ 'ui', 'gene_boxplot', 'jquery.dataTables' ]
%rebase base nav_tree=nav_tree, scripts=scripts

<link rel="stylesheet" href="/style/datatable.css" type="text/css">

<div id="main">
<h2>Taylor et al. - Integrative Genomic Profiling of Human Prostate Cancer</h2>

<h3>Genes with most often aberrated exons (primary PCa vs. BPH)</h3>
<table id="most_aberrated_primary_vs_bph" class="display"></table>

<p>
</p>

<h3>Genes with most often aberrated exons (metastasized vs. primary PCa)</h3>
<table id="most_aberrated_meta_vs_primary" class="display"></table>

</div>

<script type="text/javascript">
$(document).ready(function() {
	
	var show_table = function(id, data) {
		$.get(data, function(d) {
			var lines = d.split('\n');
			lines = lines.map(function(v) { return v.split('\t'); });
			
			for (var k = lines.length - 1; k >= 0 && lines[k] == ''; k--)
				lines.pop();
				
			for (var k = 1; k < lines.length; k++) {
				lines[k][0] = '<' + 'a href="/betastasis/prostate_cancer/taylor/alt_splicing?gene=' + lines[k][0] +
					'">' + lines[k][0] + '</a>';
			}

			
			var headers = [];
			for (var k = 0; k < lines[0].length; k++)
				headers[k] = { 'sTitle': lines[0][k] };
				
			for (var k = 0; k < lines[1].length; k++) {
				if (lines[1][k].indexOf('%') != -1)
					headers[k]['sType'] = 'percent';
			}
			
			$(id).dataTable({
				'sScrollY': '220px', 'bFilter': false,
				'bPaginate': false, 'aaData': lines.slice(1),
				'aoColumns': headers,
				'aaSorting': [[2, 'desc']]});
		});

	};
	
	show_table('#most_aberrated_primary_vs_bph',
		'/betadata/taylor_splice/most_aberrated_primary_vs_bph.txt');
	show_table('#most_aberrated_meta_vs_primary',
		'/betadata/taylor_splice/most_aberrated_metastas_vs_primary.txt');

});
</script>

