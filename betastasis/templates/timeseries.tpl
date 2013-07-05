%rebase base path=path, url=url, scripts='ui_d3,timeseries_plot'

%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Time series plot</h3>
<div id="fig"></div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Gene: <input type="text" id="GeneSelect" value="" />
</p><br>

<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}
<br><br>

<b>Export</b><br>
<button id="export_table">Table</button>
<button id="export_svg">SVG</button>
<br>
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
	var vis_data = {};
	
	var scale = 'log2';
	
	data.gene = '{{gene}}';
	var req_gene = gup('gene');
	if (req_gene.length > 0) data.gene = req_gene;
	$('#GeneSelect').val(data.gene);
	
	var vis = new TimeseriesPlot('#fig');
	
	vis.xlog = {{xlog if 'xlog' in locals() else 'false'}};
	
	var export_table = function() {
		var tabular = 'Condition\tTime\t' + data.gene + '_expr\n';
		for (var k = 0; k < data.expr.length; k++) {
			tabular += data.sample_type[k] + '\t' +
				data.timepoint[k] + '\t' + data.expr[k].toPrecision(3) + '\n';
		}
		
		download_data_uri(tabular, data.gene + '_expr_timeseries.txt',
			'application/octet-stream');
	}
	
	var downloaded = function() {
		if (!('expr' in data && 'timepoint' in data && 'sample_type' in data))
			return;
			
		vis_data = {};
		var expr = data.expr.slice(0);
		
		// Normalize expression values using the control group.
		// FIXME: Here we assume that there is only one timepoint for the
		// control samples.
		var control_expr = expr.filter(function(d, i) {
			return data.sample_type[i] == '{{control_group}}';
		});
		
		if (control_expr == undefined || control_expr.length == 0) {
			var control_mean = 0;
			vis.control_variance_per_N = 0;
		} else {
			var control_mean = d3.mean(control_expr);
			vis.control_variance_per_N = science.stats.variance(control_expr) /
				control_expr.length;
		}
		
		for (var k = 0; k < expr.length; k++) {
			expr[k] -= control_mean;
		}
		
		var timepoint = data.timepoint.map(parseFloat);

		vis.max_time = d3.max(timepoint);
		vis.max_time = vis.max_time + vis.max_time / 10;
		
		vis.min_expr = d3.min(expr);
		vis.max_expr = d3.max(expr);
		var dy = (vis.max_expr - vis.min_expr) / 10;
		vis.min_expr = vis.min_expr - dy;
		vis.max_expr = vis.max_expr + dy;
		
		
		
		// Construct a tree of the samples: first level represents conditions,
		// the second level time points.
		for (var k = 0; k < data.sample_type.length; k++) {
			if (data.sample_type[k] == '{{control_group}}')
				continue;
			
			var cond = vis_data[data.sample_type[k]];
			if (cond == undefined) {
				vis_data[data.sample_type[k]] = cond = {};
			}
			
			var time = cond[timepoint[k]];
			if (time == undefined) {
				cond[timepoint[k]] = time = [];
			}
			
			time.push(expr[k]);
		}

		vis.set_data(vis_data);
	}
	
	var download_gene = function(gene) {
		$.getJSON(data_root + '/expr/' + gene[0].toLowerCase() + '/'
			+ gene + '.json', function(d) {
			data.gene = gene; data.expr = d['data']; downloaded();
		});
	}
	
	download_gene(data.gene);
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function() {
		export_svg(data.gene + '_gene_expr_timeseries.svg'); });
	
	$.getJSON(data_root + '/clinical/Timepoint.json',
		function(d) { data.timepoint = d['data']; downloaded(); });
	$.getJSON(data_root + '/clinical/Sample_type.json',
		function(d) { data.sample_type = d['data']; downloaded(); });
	$.getJSON(data_root + '/expr/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			delete data['expr'];
			download_gene(val);
		});
	});
});
</script>

