<h3>{{'Kaplan-Meier recurrence plot' if 'recurrence' in locals() else 'Kaplan-Meier survival plot'}}</h3>
<div id="fig"></div>
</div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count"></span><br>
Censored: <span id="censored_count">-</span><br>
Gene: <input type="text" id="GeneSelect" value="" />
</p><br>

<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br><br>

<b>Statistics</b><br>
<div id="stats_div"></div>
</div>

<div id="sample_groups" class="sidebox"></div>
</div>


<script type="text/javascript">
var data_root = '{{data_root}}';

$(document).ready(function() {
	var default_gene = '{{gene if 'gene' in locals() else 'TP53'}}';
	var gene = gup('gene');
	if (gene.length == 0) { gene = default_gene; }
	$('#GeneSelect').val(gene);
	
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/groups_new.json';
	
	var vis = new GeneSurvival('fig');
	
	var data = new Object;
	var filt = new Object;
	
	vis.setStatsDiv('#stats_div');
	
	%if 'recurrence' in locals():
		vis.x_label = 'Months since treatment started';
		vis.y_label = 'Recurrence free (%)';
	%end

	var group_selector = new GroupSelector('#sample_groups', group_file,
		function() { downloaded(); });
	
	var filter_samples = function() {
		var whitelist = {};
		for (var s in group_selector.selected_samples)
			whitelist[group_selector.selected_samples[s]] = 1;
		
		var valid = function(val, idx) {
			return data.expr[idx] != 0 && data.survival[idx] != -1 &&
				data.sample_id[idx] in whitelist;
		};
		
		filt.expr = data.expr.filter(valid);
		filt.survival = data.survival.filter(valid);
		filt.censored = data.censored.filter(valid);
		
		// Update the sample count display.
		$('#sample_count').html(filt.expr.length + ' / ' + data.expr.length);
		$('#censored_count').html(filt.censored.reduce(
			function(a, b) { return b == 1 ? a + 1 : a; }) + ' / ' +
			filt.censored.length);
	}
	
	var downloaded = function() {
		if ('expr' in data && 'survival' in data && 'censored' in data &&
			'sample_id' in data && group_selector.selected_samples.length > 0) {
			filter_samples();
			vis.setData(filt.expr, filt.survival, filt.censored);
		}
	}
	
	var download_gene = function(gene) {
		$.getJSON(data_root + '/expr/' + gene[0].toLowerCase() + '/'
			+ gene + '.json', function(d) {
			data.expr = d['data']; downloaded();
		});
	}
	
	download_gene(gene);

	%if 'recurrence' in locals():
		$.getJSON(data_root + '/clinical/Recurrence-free time.json',
			function(d) { data.survival = d['data']; downloaded(); });
		$.getJSON(data_root + '/clinical/Recurrence censored.json', 
			function(d) { data.censored = d['data']; downloaded(); });
	%else:
		$.getJSON(data_root + '/clinical/Survival time.json',
			function(d) { data.survival = d['data']; downloaded(); });
		$.getJSON(data_root + '/clinical/Survival censored.json', 
			function(d) { data.censored = d['data']; downloaded(); });
	%end

	$.getJSON(data_root + '/clinical/sample_id.json', function(d) {
		data.sample_id = d['data']; downloaded();
	});
	$.getJSON(data_root + '/expr/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			delete data['expr'];
			download_gene(val);
		});
	});
});
</script>

