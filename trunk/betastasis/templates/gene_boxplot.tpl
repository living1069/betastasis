%rebase base path=path, url=url, scripts='ui,gene_boxplot'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Gene expression boxplot</h3>
<div id="fig"></div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
Gene: <input type="text" id="GeneSelect" value="" /><br/><br/>

Boxplot style:
<div id="plot_type_radio">
	<input type="radio" id="box_only" name="plot_type_radio" /><label for="box_only">Boxplot</label>
	<input type="radio" id="box_dots" name="plot_type_radio" /><label for="box_dots">Boxplot + Dots</label>
	<input type="radio" id="dots_only" name="plot_type_radio" /><label for="dots_only">Dots</label>
</div><br>

<b>Export</b><br>
<button id="export_table">Table</button>
<button id="export_svg">SVG</button><br>
<br>

<b>Platforms</b><br>
{{!platform if 'platform' in locals() else '-'}}<br>
</p>
</div>

</div>

<script type="text/javascript">
var data_root = '{{data}}';

$(document).ready(function() {
	
	var data = {gene: ''};
	var filt = [];
	
	var default_gene = gup('gene');
	if (default_gene.length == 0)
		default_gene = '{{gene if 'gene' in locals() else 'TP53'}}';
	$('#GeneSelect').val(default_gene);
	
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/box_groups.json';

	var data_downloaded = function() {
		if (!('sample_ids' in data && 'groups' in data && 'expr' in data))
			return;
			
		filt.length = 0;
		for (var g in data.groups) {
			var gsamples = data.groups[g];
			
			var whitelist = {};
			for (var s in gsamples)
				whitelist[gsamples[s]] = 1;
				
			var validate = function(val, idx) {
				return data.expr[idx] != 0 && data.sample_ids[idx] in whitelist;
			}
				
			var gexpr = data.expr.filter(validate);
			var filt_sample_ids = data.sample_ids.filter(validate);
			
			var quart = bmath.quartiles(gexpr);
				
			filt.push({
				subtype: g,
				min: pv.min(gexpr),
				max: pv.max(gexpr),
				median: quart.median,
				lq: quart.lq,
				uq: quart.uq,
				expr_values: gexpr,
				sample_ids: filt_sample_ids
			});
		}
		vis.gene_name = data.gene;
		vis.set_data(filt);
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
	
	var export_table = function() {
		var tabular = '';
		for (var g = 0; g < filt.length-1; g++)
			tabular += filt[g].subtype + '\t';
		tabular += filt[filt.length-1].subtype + '\n';
		
		for (var k = 0; ; k++) {
			var done = true;
			for (var g = 0; g < filt.length; g++) {
				if (filt[g].expr_values.length <= k) continue;
				tabular += filt[g].expr_values[k].toFixed(2);
				if (g < filt.length - 1) tabular += '\t';
				done = false;
			}
			tabular += '\n';
			if (done) break;
		}
		
		download_data_uri(tabular, data.gene + '_expr_boxplot.txt',
			'application/octet-stream');
	}
	
	$('#plot_type_radio').buttonset();
	$('#box_only').attr('checked', true);
	$('#box_only').button('refresh');

	var vis = new ExprBoxplot('fig');

	$("#box_only").button().click(function() {
		vis.set_box_visibility(true);
		vis.set_dots_visibility(false);
	});
	$("#box_dots").button().click(function() {
		vis.set_box_visibility(true);
		vis.set_dots_visibility(true);
	});
	$("#dots_only").button().click(function() {
		vis.set_box_visibility(false);
		vis.set_dots_visibility(true);
	});
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function() {
		export_svg(data.gene + '_expr_boxplot.svg') });

	$.getJSON(group_file, function(d) {
		data.groups = d; data_downloaded();
	});

    load_gene_data(default_gene);
	
    $.getJSON(data_root + '/expr/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			delete data['expr']; load_gene_data(val);
		});
	});
	
	$.getJSON(data_root + '/clinical.json', function(d) {
		data.sample_ids = d['sample_id']; data_downloaded();
	});

});
</script>

