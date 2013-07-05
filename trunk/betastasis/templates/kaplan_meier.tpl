%rebase base path=path, url=url, scripts='ui,gene_survival'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Kaplan-Meier plot</h3>
<div id="fig"></div>
</div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count"></span><br>
Censored: <span id="censored_count">-</span><br>
{{feature_type if 'feature_type' in locals() else 'Gene'}}: <input type="text" id="GeneSelect" value="" />
</p><br>

<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br><br>

<b>Preset thresholds</b><br>
<button id="threshold_25th">25%</button>
<button id="threshold_median">Median</button>
<button id="threshold_75th">75%</button>
<br><br>

<b>Statistics</b><br>
<div id="stats_div"></div><br>

<b>Export</b><br>
<button id="export_table">Table</button>
<button id="export_svg">SVG</button>
<br>
</div>


<div id="sample_groups" class="sidebox"></div>
</div>


<script type="text/javascript">
var data_root = '{{data}}';

$(document).ready(function() {
	
	var data = {};
	var filt = {};
	
	var group_file = '{{groups}}';
	
	var vis = new GeneSurvival('fig');
	vis.setStatsDiv('#stats_div');
	
	var group_selector = new GroupSelector('#sample_groups', group_file,
		function() { downloaded(); });
	
	var export_table = function() {
		var tabular = 'Sample\tExpression\tSurvival_time\tCensored\n';
		for (var k = 0; k < filt.sample_id.length; k++) {
			tabular += filt.sample_id[k] + '\t' + filt.expr[k] + '\t' +
				filt.survival[k] + '\t' + filt.censored[k] + '\n';
		}
		
		download_data_uri(tabular, data.feature + '_kaplan_meier.txt', 
			'application/octet-stream');
	}
	
	var filter_samples = function() {
		var whitelist = {};
		
		for (var s in group_selector.selected_samples)
			whitelist[group_selector.selected_samples[s]] = 1;
		
		var valid = function(val, idx) {
			return data.expr[idx] != 0 && data.survival[idx] != -1 &&
				data.sample_id[idx] in whitelist;
		};
		
		filt.sample_id = data.sample_id.filter(valid);
		filt.expr = data.expr.filter(valid);
		filt.survival = data.survival.filter(valid);
		filt.censored = data.censored.filter(valid);
		
		// Update the sample count display.
		$('#sample_count').html(filt.expr.length + ' / ' + data.expr.length);
		$('#censored_count').html(filt.censored.filter(
			function(a) { return a == 1; }).length + ' / ' +
			filt.censored.length);
	}
	
	var downloaded = function() {
		if ('expr' in data && 'survival' in data && 'censored' in data &&
			'sample_id' in data && group_selector.selected_samples.length > 0) {
			filter_samples();
			vis.setData(filt.expr, filt.survival, filt.censored);
		}
	}
	
	var download_feature = function(feature) {
		// Without this check, clicking on "Median" after changing the displayed
		// gene yields bad behavior.
		if (data.feature == feature)
			return;
		
		data.feature = feature;
		$.getJSON(data_root + '/expr/' + feature[0].toLowerCase() + '/'
			+ feature + '.json', function(d) {
			data.expr = d['data']; downloaded();
		});
	}
	
	var feature = gup('feature');
	if (feature.length == 0) { feature = '{{feature}}'; }
	$('#GeneSelect').val(feature);

	download_feature(feature);
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function() {
		export_svg(data.feature + '_kaplan_meier.svg'); });
	
	$('#threshold_25th').button().click(function() {
		vis.setExprThreshold(bmath.quartiles(filt.expr).lq);
	});
	$('#threshold_median').button().click(function() {
		vis.setExprThreshold(pv.median(filt.expr));
	});
	$('#threshold_75th').button().click(function() {
		vis.setExprThreshold(bmath.quartiles(filt.expr).uq);
	});
	
	$.getJSON(data_root + '/clinical.json', function(d) {
		if ('survival_time' in d) {
			data.survival = d['survival_time'];
			data.censored = d['survival_censored'];
		} else if ('progression_time' in d) {
			data.survival = d['progression_time'];
			data.censored = d['progression_censored'];
			
			vis.x_label = 'Months since treatment started';
			vis.y_label = 'Recurrence free (%)';
		}
		data.sample_id = d['sample_id'];
		downloaded();
	});

	$.getJSON(data_root + '/expr/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			delete data['expr'];
			download_feature(val);
		});
	});
});
</script>

