%rebase base path=path, url=url, scripts='ui,barplot'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>{{feature_type if 'feature_type' in locals() else 'Gene'}} expression bar plot</h3>
<div id="fig"></div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
{{feature_type if 'feature_type' in locals() else 'Gene'}}: <input type="text" id="GeneSelect" value="" /><br><br>

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

<div id="sample_groups" class="sidebox"></div>

</div>


%if defined('description'):
<h3>Description</h3>
{{description}}
%end
</div>

	  
<script type="text/javascript">
var data_root = '{{data}}';

$(document).ready(function() {
	
	var data = {gene: ''};
	var filtered;
	
	var default_gene = gup('gene');
	if (default_gene.length == 0)
		default_gene = '{{gene if 'gene' in locals() else 'TP53'}}';
	$('#GeneSelect').val(default_gene);
	
	// Initialize a group selector only if the user has specified a JSON
	// file containing sample groupings.
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	var group_selector = undefined;
	if (group_file.length > 0) {
		group_selector = new GroupSelector('#sample_groups', group_file,
			function() { data_downloaded(); });
	}
	
	$('#scale_radio').buttonset();
	$('#natural_scale').attr('checked', true);
	$('#natural_scale').button('refresh');
	
	var data_downloaded = function() {
		if (!('sample_ids' in data && 'expr' in data))
			return;
		
		// Filter samples only if a group JSON file has been provided.
		if (typeof group_selector !== 'undefined') {
			if (group_selector.ready == false) return;
			
			var valid = function(val, idx) {
				return group_selector.selected_samples
					.indexOf(data.sample_ids[idx]) > -1;
			};
		
			filtered = {};
			filtered.sample_ids = data.sample_ids.filter(valid);
			filtered.expr = data.expr.filter(valid);
		} else {
			filtered = data;
		}

		vis.gene_name = data.gene;
		vis.set_data(filtered.expr, filtered.sample_ids);
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
		var tabular = 'Sample\t' + data.gene + ' ' +
			(vis.scale == 'log-2' ? 'log-2' : '') + ' expr\n';
		for (var k = 0; k < filtered.sample_ids.length; k++) {
			var val = filtered.expr[k];
			if (vis.scale == 'natural') val = Math.pow(2, val);
			tabular += filtered.sample_ids[k] + '\t' + val.toFixed(3) + '\n';
		}
		
		download_data_uri(tabular, data.gene + '_expr.txt',
			'application/octet-stream');
	}

	var vis = new Barplot('fig');
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function() {
		export_svg(data.gene + '_expr_barplot.svg'); });
	
	$("#natural_scale").button().click(function() {
		vis.scale = 'natural';
		vis.set_data(filtered.expr, filtered.sample_ids);
	});
	$("#log2_scale").button().click(function() {
		vis.scale = 'log-2';
		vis.set_data(filtered.expr, filtered.sample_ids);
	});

    load_gene_data(default_gene);
	
    $.getJSON(data_root + '/expr/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			load_gene_data(val);
		});
	});
	
	$.getJSON(data_root + '/clinical.json', function(d) {
		data.sample_ids = d['sample_id']; data_downloaded();
	});

});
</script>

