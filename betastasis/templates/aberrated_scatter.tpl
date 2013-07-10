%rebase base path=path, url=url, scripts='ui_d3,scatterplot.d3,fisheye'

%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Most aberrated genes</h3>
<div id="fig"></div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
</p><br>

<b>Platforms</b><br>
{{!platform if 'platform' in locals() else '-'}}
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
	
	var dots = {};
	
	var vis = new Scatterplot('#fig');
	
	var export_table = function() {
		/*
		var tabular = 'Condition\tTime\t' + data.gene + '_expr\n';
		for (var k = 0; k < data.expr.length; k++) {
			tabular += data.sample_type[k] + '\t' +
				data.timepoint[k] + '\t' + data.expr[k].toPrecision(3) + '\n';
		}
		
		download_data_uri(tabular, data.gene + '_expr_timeseries.txt',
			'application/octet-stream');
		*/
	}
	
	var size_scale = d3.scale.linear()
		.domain([-3,0,3]).range([6,1,6]).clamp(true);
	var saturation_scale = d3.scale.linear()
		.domain([-3,0,3]).range([100,0,100]).clamp(true);
	var light_scale = d3.scale.linear()
		.domain([-3,0,3]).range([50,95,50]).clamp(true);
	
	var downloaded = function(d) {
		dots.x = d.values[1];
		dots.y = d.values[2];
		dots.radius = d.values[0].map(
			function(v) { return size_scale(v) + 'pt'; });
		dots.color = d.values[0].map(function(v) {
			return 'hsl(' + (v < 0 ? 240 : 0) + ',' +
				saturation_scale(v) + '%,' + light_scale(v) + '%)';
		});
			
		dots.label = d.labels;
		vis.set_data(dots);
	}
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function() {
		export_svg('aberration_scatter.svg'); });
	
	d3.json(data_root, downloaded);
	
});
</script>

