%rebase base path=path, url=url, scripts='ui,barplot'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Exon expression bar plot</h3>
<div id="fig"></div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
Gene: <input type="text" id="gene_select" value="" /><br>
Exon: <span id="exon_select"></span><br><br>

Expression scale:
<div id="scale_radio">
	<input type="radio" id="natural_scale" name="scale_radio" /><label for="natural_scale">Natural</label>
	<input type="radio" id="log2_scale" name="scale_radio" /><label for="log2_scale">Log-2</label>
</div><br>

<b>Export</b><br>
<button id="export_table">Table</button>
<button id="export_svg">SVG</button><br>
<br>



<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br>
</p>

</div>
</div>

%if defined('description'):
<h3>Description</h3>
{{description}}
%end
	  
<script type="text/javascript">
var data_root = '{{data}}';

$(document).ready(function() {
	
	var data = {};
	
	var req_gene = gup('gene');
	data.gene = (req_gene.length > 0) ? req_gene : '{{gene}}';
	
	var req_exon = gup('exon');
	data.default_exon = (req_exon.length > 0) ? req_exon : '{{exon}}';
		
	$('#gene_select').val(data.gene);
	
	$('#scale_radio').buttonset();
	$('#natural_scale').attr('checked', true);
	$('#natural_scale').button('refresh');
	
	var data_downloaded = function() {
		if (!('sample_ids' in data && 'expr' in data))
			return;
			
		var exon_expr = data.expr[data.exon_idx].slice(0);
		
		vis.gene_name = data.gene;
		vis.set_data(exon_expr, data.sample_ids);
	}

	var download_exon_data = function(gene, default_exon) {
		$.getJSON(data_root + '/exon_expr/' + gene[0].toLowerCase() + '/'
			+ gene + '.json', function(d) {
			
			// Update the exon list.
			data.exons = d['exons'];
			data.exon_idx = 0;
			for (var k in data.exons) {
				if (default_exon == data.exons[k])
					data.exon_idx = k;
			}
			
			// Grab the expression data and update the visualization.
			data.expr = d['data'];
			data.gene = gene;
			
			update_exon_list();
			
			data_downloaded();
		});
	}
	
	var update_exon_list = function() {
		var str = '';
		var exons = data.exons;
		for (var e = 0; e < exons.length; e++) {
			str += '<a id="exon_' + e + '" href="">';
			str += (e == data.exon_idx) ? 
				('<b>' + exons[e] + '</b>') : exons[e];
			str += '</a>';
			if (e < exons.length-1) str += ' | ';
		}
		$('#exon_select').html(str);
		
		for (var e = 0; e < exons.length; e++) {
			$('#exon_' + e).click((function(ex) {
				return function() {
					if (ex == data.exon_idx) return false;
					data.exon_idx = ex;
					data_downloaded();
					update_exon_list();
					return false;
				}
			})(e));
		}
	}
	
	var export_table = function() {
		var tabular = 'Sample';
		
		for (var k = 0; k < data.exons.length; k++)
			tabular += '\t' + data.gene + '[' + data.exons[k] + ']';
		tabular += '\n';
		
		for (var s = 0; s < data.sample_ids.length; s++) {
			tabular += data.sample_ids[s];
			for (var k = 0; k < data.exons.length; k++) {
				var val = data.expr[k][s];
				if (vis.scale == 'natural')
					val = Math.pow(2, val);
				tabular += '\t' + val.toFixed(3);
			}
				
			tabular += '\n';
		}
		
		download_data_uri(tabular, data.gene + '_exon_expr.txt', 
			'application/octet-stream');
	}

	var vis = new Barplot('fig');
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function() {
		export_svg(data.gene + '_exon_expr_barplot.svg'); });
	
	$("#natural_scale").button().click(function() {
		vis.scale = 'natural'; data_downloaded(); });
	$("#log2_scale").button().click(function() {
		vis.scale = 'log-2'; data_downloaded(); });

    download_exon_data(data.gene, data.default_exon);
	
    $.getJSON(data_root + '/expr/features.json', function(d) {
		$('#gene_select').omnicomplete(d['features'], function(val) {
			download_exon_data(val);
		});
	});
	
	$.getJSON(data_root + '/clinical/sample_id.json', function(d) {
		data.sample_ids = d['data']; data_downloaded();
	});

});
</script>

