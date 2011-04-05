
<h3>Gene expression boxplot</h3>
<div id="fig"></div>
</div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
Gene: <input type="text" id="GeneSelect" value="" /><br>
Exon: <span id="exon_select"></span><br><br>

Boxplot style:
<div id="plot_type_radio">
	<input type="radio" id="box_only" name="plot_type_radio" /><label for="box_only">Boxplot</label>
	<input type="radio" id="box_dots" name="plot_type_radio" /><label for="box_dots">Boxplot + Dots</label>
	<input type="radio" id="dots_only" name="plot_type_radio" /><label for="dots_only">Dots</label>
</div><br>

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
	
	var default_exon = gup('exon');
	
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/box_groups.json';
	
	var data_downloaded = function() {
		if (!('sample_ids' in data && 'groups' in data && 'expr' in data))
			return;
			
		var e = data.selected_exon;
		
		var filt = [];
		for (var g in data.groups) {
			var gsamples = data.groups[g];
			
			var whitelist = {};
			for (var s in gsamples)
				whitelist[gsamples[s]] = 1;
				
			var validate = function(val, idx) {
				return data.expr[e][idx] != 0 &&
					data.sample_ids[idx] in whitelist;
			}
				
			var gexpr = data.expr[e].filter(validate);
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

	var update_exon_list = function() {
		var str = '';
		var exons = data.exons;
		for (var e = 0; e < exons.length; e++) {
			str += '<a id="exon_' + e + '" href="">';
			str += (e == data.selected_exon) ? 
				('<b>' + exons[e] + '</b>') : exons[e];
			str += '</a>';
			if (e < exons.length-1) str += ' | ';
		}
		$('#exon_select').html(str);
		
		for (var e = 0; e < exons.length; e++) {
			$('#exon_' + e).click((function(ex) {
				return function() {
					if (ex == data.selected_exon) return false;
					data.selected_exon = ex;
					data_downloaded();
					update_exon_list();
					return false;
				}
			})(e));
		}
	}

	var load_gene_data = function(gene_name) {
		if (gene_name == data.gene) return;
		$.getJSON(data_root + '/exon_splice/' + gene_name[0].toLowerCase() + '/'
			+ gene_name + '.json', function(d) {
				
			// Update the exon list.
			data.exons = d['exons'];
			data.selected_exon = 0;
			for (var k in data.exons) {
				if (default_exon == data.exons[k])
					data.selected_exon = k;
			}
			default_exon = '';   // Only use the default exon once.
			update_exon_list();
			
			// Grab the expression data and update the visualization.
			data.expr = d['data'];
			data.gene = gene_name;

			data_downloaded();
		});
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

	$.getJSON(group_file, function(d) {
		data.groups = d; data_downloaded();
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

