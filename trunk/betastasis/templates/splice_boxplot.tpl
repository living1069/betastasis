%rebase base path=path, url=url, scripts='ui,gene_boxplot'

%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Splicing index boxplot</h3>
<div id="fig"></div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
Gene: <input type="text" id="GeneSelect" value="" /><br><br>

Boxplot style:
<div id="plot_type_radio">
	<input type="radio" id="box_only" name="plot_type_radio" /><label for="box_only">Boxplot</label>
	<input type="radio" id="box_dots" name="plot_type_radio" /><label for="box_dots">Boxplot + Dots</label>
	<input type="radio" id="dots_only" name="plot_type_radio" /><label for="dots_only">Dots</label>
</div><br>

<b>Platforms</b><br>
{{!platform if 'platform' in locals() else '-'}}
</p>
</div>

<div id="test_groups" class="sidebox"></div>
<div id="ref_groups" class="sidebox"></div>

</div>


	  
<script type="text/javascript">
var data_root = '{{data}}';

$(document).ready(function() {
	
	var data = {gene: ''};
	var test_filt = new Object;
	var ref_filt = new Object;
	
	var default_gene = gup('gene');
	if (default_gene.length == 0)
		default_gene = '{{gene if 'gene' in locals() else 'TP53'}}';
	$('#GeneSelect').val(default_gene);
	
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/groups_new.json';
	
	var test_group_selector = new GroupSelector('#test_groups', group_file,
		function() { data_downloaded(); });
	test_group_selector.title = 'Test group';
	test_group_selector.default_group = '{{default_test_group}}';
	
	var ref_group_selector = new GroupSelector('#ref_groups', group_file,
		function() { data_downloaded(); });
	ref_group_selector.title = 'Reference group';
	ref_group_selector.default_group = '{{default_ref_group}}';
	
	var data_downloaded = function() {
		if (!('sample_ids' in data && 'expr' in data && 
			test_group_selector.selected_samples.length > 0 &&
			ref_group_selector.selected_samples.length > 0))
			return;
			
		// Pick test group samples.
		var whitelist = {};
		for (var s in test_group_selector.selected_samples)
			whitelist[test_group_selector.selected_samples[s]] = 1;
		
		var valid = function(val, idx) {
			return data.sample_ids[idx] in whitelist;
		};
		
		test_filt.expr = [];
		for (var e = 0; e < data.expr.length; e++)
			test_filt.expr[e] = data.expr[e].filter(valid);
		test_filt.sample_ids = data.sample_ids.filter(valid);
		
		// Pick reference group samples.
		var whitelist = {};
		for (var s in ref_group_selector.selected_samples)
			whitelist[ref_group_selector.selected_samples[s]] = 1;
		
		var valid = function(val, idx) {
			return data.sample_ids[idx] in whitelist;
		};
		
		ref_filt.expr = [];
		for (var e = 0; e < data.expr.length; e++)
			ref_filt.expr[e] = data.expr[e].filter(valid);
			
		var groups = [];
		
		// Calculate differential splicing between the two groups.
		for (var e = 0; e < ref_filt.expr.length; e++) {
			var pad = 100;
			var ref_median = pv.median(ref_filt.expr[e], function(d) {
				if (d == 0) { pad = -pad; return pad; } else { return d; }
			});
			
			var gexpr = [];
			var gsamples = [];
			
			var exon_expr = test_filt.expr[e];
			for (var s = 0; s < exon_expr.length; s++) {
				if (exon_expr[s] == 0) continue;
				gexpr.push(exon_expr[s] - ref_median);
				gsamples.push(test_filt.sample_ids[s]);
			}
			
			var quart = bmath.quartiles(gexpr);
			
			groups.push({
				subtype: data.exons[e],
				min: pv.min(gexpr),
				max: pv.max(gexpr),
				median: quart.median,
				lq: quart.lq,
				uq: quart.uq,
				expr_values: gexpr,
				sample_ids: gsamples
			});
		}
		vis.gene_name = data.gene;
		vis.set_data(groups);
	}

	var load_gene_data = function(gene_name) {
		if (gene_name == data.gene) return;
		$.getJSON(data_root + '/exon_splice/' + gene_name[0].toLowerCase() + '/'
			+ gene_name + '.json', function(d) {
			data.exons = d['exons'];
			data.expr = d['data'];
			data.gene = gene_name;
			data_downloaded();
		});
	}
	
	$('#plot_type_radio').buttonset();
	$('#box_only').attr('checked', true);
	$('#box_only').button('refresh');

	var vis = new ExprBoxplot('fig');
	vis.diag_labels = false;

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

    load_gene_data(default_gene);
	
    $.getJSON(data_root + '/exon_splice/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			load_gene_data(val);
		});
	});
	
	$.getJSON(data_root + '/clinical/sample_id.json', function(d) {
		data.sample_ids = d['data']; data_downloaded();
	});

});
</script>

<h3>Description</h3>
<p>This visualization displays relative exon expression levels between two groups of samples. Exons are numbered along the X-axis of the boxplot, and log-2 differential expression is shown on the Y-axis. For each exon, a median expression level is calculated based on the reference sample group. Then the expression of each exon in each sample in the test group is compared against the reference expression level.</p>

</div>

