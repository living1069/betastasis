%rebase base path=path, url=url, scripts='ui,scatterplot'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Gene expression scatterplot</h3>
<div id="fig"></div>
</div>
	  
<div id="sidebar">

<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count"></span><br>
Gene X: <input type="text" id="gene_x_select" value="" /><br/>
Gene Y: <input type="text" id="gene_y_select" value="" /><br/><br/>

<div id="plot_scale">
	<input type="radio" id="natural_scale" name="plot_scale" /><label for="natural_scale">Natural</label>
	<input type="radio" id="log2_scale" name="plot_scale" /><label for="log2_scale">Log-2</label>
</div><br>

<b>Statistics</b><br>
Pearson = <span id="pearson"></span><br><br>

<b>Curve fitting</b>
<div id="fit_stats_div"></div><br>

<b>Export</b><br>
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

$(document).ready(function() {
	
	var data = {gene_x: '', gene_y: ''};
	var filt = {};
	
	data.gene_x = gup('gene_x');
	if (data.gene_x.length == 0)
		data.gene_x = '{{gene_x}}';
	$('#gene_x_select').val(data.gene_x);
	
	data.gene_y = gup('gene_y');
	if (data.gene_y.length == 0)
		data.gene_y = '{{gene_y}}';
	$('#gene_y_select').val(data.gene_y);
	
	var group_file = '{{groups}}';
	var group_selector = new GroupSelector('#sample_groups', group_file,
		function() { data_downloaded(); });
	
	var vis = new Scatterplot('fig');
	
	var data_downloaded = function() {
		if (!('sample_ids' in data && 'expr_x' in data && 'expr_y' in data &&
			group_selector.selected_samples.length > 0))
			return;
			
		// Only pick samples that belong to the selected sapmle group.
		// Also filter out any samples with missing values for the gene.
		var whitelist = {};
		for (var s in group_selector.selected_samples)
			whitelist[group_selector.selected_samples[s]] = 1;
		
		var valid = function(val, idx) {
			return data.expr_x[idx] != 0 && data.expr_y[idx] != 0 &&
				data.sample_ids[idx] in whitelist;
		};
		
		filt.xval = data.expr_x.filter(valid);
		filt.yval = data.expr_y.filter(valid);
		filt.sample_ids = data.sample_ids.filter(valid);
		
		// Update the sample count display.
		$('#sample_count').html(filt.xval.length + ' / ' +
			data.expr_x.length);

		vis.xlabel = data.gene_x;
		vis.ylabel = data.gene_y;
		vis.set_data(filt);
		
		$('#pearson').html(vis.pearson.toPrecision(3));
		
		$('#fit_stats_div').html(
			'<button id="fit_button">Least squares fit</button><br>');
		$('#fit_button').button();
			
		$('#fit_button').click(function() {
			coeff = vis.fit_least_squares();
			$('#fit_stats_div').html(
				'Intercept = ' + coeff[0].toPrecision(5) + '<br>' +
				'1st order = ' + coeff[1].toPrecision(5));
		});
	}

	var load_gene_data = function() {
		delete data['expr_x']; 
		delete data['expr_y'];
		
		$.getJSON(data_root + '/expr/' + data.gene_x[0].toLowerCase() + '/'
			+ data.gene_x + '.json', function(d) {
			data.expr_x = d['data'];
			data_downloaded();
		});
		
		$.getJSON(data_root + '/expr/' + data.gene_y[0].toLowerCase() + '/'
			+ data.gene_y + '.json', function(d) {
			data.expr_y = d['data'];
			data_downloaded();
		});
	}

	$('#plot_scale').buttonset();
	$('#log2_scale').attr('checked', true);
	$('#log2_scale').button('refresh');
	
	$('#export_svg').button().click(function() {
		export_svg(data.gene_y + '_vs_' + data.gene_x + '_scatter.svg'); });

	$("#natural_scale").button().click(function() {
		vis.scale = 'natural'; data_downloaded();
	});
	$("#log2_scale").button().click(function() {
		vis.scale = 'log2'; data_downloaded();
	});

	$.getJSON(data_root + '/clinical.json', function(d) {
		data.sample_ids = d['sample_id']; data_downloaded();
	});

    $.getJSON(data_root + '/expr/features.json', function(d) {
		$('#gene_x_select').omnicomplete(d['features'], function(val) {
			if (data.gene_x != val) {
				data.gene_x = val; load_gene_data();
			}
		});
		
		$('#gene_y_select').omnicomplete(d['features'], function(val) {
			if (data.gene_y != val) {
				data.gene_y = val; load_gene_data();
			}
		});
	});
	
	load_gene_data();
});
</script>

