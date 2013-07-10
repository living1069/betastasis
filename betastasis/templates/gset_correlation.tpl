%rebase base path=path, url=url, scripts='ui,vertical_barplot'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Gene set correlation</h3>
<div id="fig"></div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count"></span><br>
Gene: <input type="text" id="ref_gene_select" value="" /><br>
Targets: <input type="text" id="target_genes_select" value="" /><br><br>

<b>Correlation type</b><br>
<div id="correlation_radio">
	<input type="radio" id="pearson_corr" name="correlation_radio" /><label for="pearson_corr">Pearson</label>
	<!--<input type="radio" id="spearman_corr" name="correlation_radio" /><label for="spearman_corr">Spearman</label>-->
</div><br>

<b>Sorting</b><br>
<div id="sorting_radio">
	<input type="radio" id="no_sorting" name="sorting_radio" /><label for="no_sorting">None</label>
	<input type="radio" id="descend_sorting" name="sorting_radio" /><label for="descend_sorting">Descend</label>
	<input type="radio" id="ascend_sorting" name="sorting_radio" /><label for="ascend_sorting">Ascend</label>
</div><br>

<b>Export</b><br>
<button id="export_table">Table</button>
<button id="export_svg">SVG</button><br>
<br>


<b>Platforms</b><br>
{{!platform if 'platform' in locals() else '-'}}<br>
</p>
</div>

<div id="sample_groups" class="sidebox"></div>
</div>
	  
<script type="text/javascript">
var data_root = '{{data}}';
var scatter_url = '{{scatter_url}}';

var pearson_corr = function(a, b) {
	var a_mean = pv.mean(a), a_stdev = pv.deviation(a);
	var b_mean = pv.mean(b), b_stdev = pv.deviation(b);
	
	var pearson = 0;
	for (var k = 0; k < a.length; k++)
		pearson += (a[k] - a_mean) * (b[k] - b_mean);
	pearson /= a_stdev * b_stdev * a.length;
	return pearson;
};

var spearman_corr = function(a, b) {
	// FIXME: Handle tied values properly.
	var aidx = a.map(function(d, idx) { return [idx, d]; });
	var bidx = b.map(function(d, idx) { return [idx, d]; });
	
	aidx = aidx.sort(function(a, b) { return a[1] - b[1]; })
	bidx = bidx.sort(function(a, b) { return a[1] - b[1]; })
	
	var arank = aidx.map(function(d) { return d[0]; });
	var brank = bidx.map(function(d) { return d[0]; });
	return pearson_corr(arank, brank);
}

$(document).ready(function() {
	
	var data = {ref_expr: [], ref_gene: '', target_expr: [], target_genes: [],
		target_gene_str: ''};
	var vis_data = {};
	
	var sorting = 'none';
	var corr_type = 'pearson';
	
	data.ref_gene = gup('ref_gene');
	if (data.ref_gene.length == 0)
		data.ref_gene = '{{ref_gene}}';
	$('#ref_gene_select').val(data.ref_gene);
	
	data.target_gene_str = gup('target_genes');
	if (data.target_gene_str.length == 0)
		data.target_gene_str = '{{target_genes}}';
	$('#target_genes_select').val(data.target_gene_str);
	
	var group_file = '{{groups}}';
	var group_selector = new GroupSelector('#sample_groups', group_file,
		function() { data_downloaded(); });
	
	var vis = new VerticalBarplot('fig');
	vis.limits = [-1, 1];
	vis.sorting = 'none';
	
	$('#correlation_radio').buttonset();
	$('#pearson_corr').attr('checked', true);
	$('#pearson_corr').button('refresh');
	
	$('#sorting_radio').buttonset();
	$('#no_sorting').attr('checked', true);
	$('#no_sorting').button('refresh');
	
	var export_table = function() {
		var tabular = 'Gene\t' + corr_type[0].toUpperCase() + 
			corr_type.substr(1) + '_with_' + data.ref_gene + '\n';
		for (var k = 0; k < vis_data.values.length; k++) {
			tabular += vis_data.labels[k] + '\t' +
				vis_data.values[k].toPrecision(5) + '\n';
		}
		
		download_data_uri(tabular, data.ref_gene + '_gset_correlation.txt',
			'application/octet-stream');
	}
	
	var data_downloaded = function() {
		if (!('sample_ids' in data && 'ref_expr' in data &&
			group_selector.selected_samples.length > 0))
			return;
			
		var labels = [];
		var values = [];
		var ref = [];
		
		for (var k = 0; k < data.target_genes.length; k++) {
			if (data.target_expr[k].length > 0) {
				labels.push(data.target_genes[k]);
				values.push(data.target_expr[k]);
			}
		}
		
		if (labels.length == 0) return;
		
		// Discard samples which do not have an expression value for one of
		// the features.
		var valid = pv.repeat([1], values[0].length);
		for (var d = 0; d < values.length; d++) {
			var ds = values[d];
			if (data.ref_expr[s] == 0) valid[s] = 0;
			for (var s = 0; s < ds.length; s++) {
				if (ds[s] == 0) valid[s] = 0;
			}
		}
		
		var samples_before_filtering = values[0].length;
		
		var whitelist = {};
		for (var s in group_selector.selected_samples)
			whitelist[group_selector.selected_samples[s]] = 1;

		var validator = function(val, idx) {
			return valid[idx] && data.sample_ids[idx] in whitelist;
		};

		for (var d = 0; d < values.length; d++)
			values[d] = values[d].filter(validator);
		ref = data.ref_expr.filter(validator);
		
		// Now we must calculate Pearson correlations between the reference
		// gene and target genes.
		var corr = [];
		for (var d = 0; d < values.length; d++) {
			if (corr_type == 'pearson') {
				corr[d] = pearson_corr(ref, values[d]);
			} else {
				corr[d] = spearman_corr(ref, values[d]);
			}
		}
		
		// Sort the data if requested by the user.
		if (sorting != 'none') {
			var to_sort = [labels, corr];
			pv.transpose(to_sort);
			
			if (sorting == 'ascend') {
				to_sort.sort(function(a, b) { return a[1] - b[1]; });
			} else if (sorting == 'descend') {
				to_sort.sort(function(a, b) { return b[1] - a[1]; });
			} else {
				alert('Invalid sort order')
			}
			
			pv.transpose(to_sort);
			labels = to_sort[0];
			corr = to_sort[1];
		}
		
		// Update the sample count display.
		$('#sample_count').html(values[0].length + ' / ' + 
			samples_before_filtering);
		
		vis.xlabel = 'Correlation with ' + data.ref_gene;

		vis_data = {labels: labels, values: corr};
		
		// Construct links to scatterplots.
		var urls = [];
		for (var k = 0; k < labels.length; k++) {
			urls[k] = scatter_url + '?gene_x=' +
				data.ref_gene + '&gene_y=' + labels[k];
		}
		vis_data['urls'] = urls;
		
		vis.set_data(vis_data);
	}

	var load_gene_data = function() {
		
		// Load the reference gene data.
		$.getJSON(data_root + '/expr/' + data.ref_gene[0].toLowerCase() + '/'
			+ data.ref_gene.toUpperCase() + '.json', function(d) {
			data.ref_expr = d['data'];
			data_downloaded();
		});
		
		// Load data for the target genes.
		var targets = data.target_gene_str.match(/[\w-]+/g);
		for (var k = 0; k < targets.length; k++)
			data.target_expr[k] = [];
			
		data.target_genes = targets;

		for (var k = 0; k < targets.length; k++) {
			var target = targets[k];
			
			$.getJSON(data_root + '/expr/' + target[0].toLowerCase() + '/'
				+ target.toUpperCase() + '.json', 
				(function(idx, feature) { return function(d) {
					data.target_expr[idx] = d['data'];
					data_downloaded();
				};})(k, target)
			);
		}
	}

	$("#pearson_corr").button().click(function() {
		corr_type = 'pearson'; data_downloaded(); });
	/*$("#spearman_corr").button().click(function() {
		corr_type = 'spearman'; data_downloaded(); });*/
	
	$('#no_sorting').button().click(function() {
		sorting = 'none'; data_downloaded(); });
	$('#descend_sorting').button().click(function() {
		sorting = 'descend'; data_downloaded(); });
	$('#ascend_sorting').button().click(function() {
		sorting = 'ascend'; data_downloaded(); });
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function() {
		export_svg(data.ref_gene + '_gset_correlation.svg'); });

    load_gene_data();
	
	$('#target_genes_select').change(function() {
		data.target_gene_str = $(this).val();
		load_gene_data();
	});
	
    $.getJSON(data_root + '/expr/features.json', function(d) {
		$('#ref_gene_select').omnicomplete(d['features'], function(val) {
			if (data.ref_gene != val) {
				data.ref_gene = val;
				load_gene_data();
			}
		});
	});
	
	$.getJSON(data_root + '/clinical/sample_id.json', function(d) {
		data.sample_ids = d['data']; data_downloaded();
	});

});
</script>

