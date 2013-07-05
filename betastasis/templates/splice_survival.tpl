%rebase base path=path, url=url, scripts='ui,gene_survival'

%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>{{'Alternative splicing Kaplan-Meier recurrence plot' if 'recurrence' in locals() else 'Alternative splicing Kaplan-Meier survival plot'}}</h3>
<div id="fig"></div>
</div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count"></span><br>
Censored: <span id="censored_count">-</span><br>
Gene: <input type="text" id="GeneSelect" value="" /><br>
Exon: <span id="exon_select"></span>
</p><br>

<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br><br>

<b>Statistics</b><br>
<div id="stats_div"></div>
</div>

<div id="test_groups" class="sidebox"></div>

</div>



<script type="text/javascript">
var data_root = '{{data}}';

$(document).ready(function() {
	var default_gene = '{{gene if 'gene' in locals() else 'TP53'}}';
	var gene = gup('gene');
	if (gene.length == 0) { gene = default_gene; }
	$('#GeneSelect').val(gene);
	
	var default_exon = gup('exon');
	
	var data = new Object;
	var test_filt = new Object;
	
	var vis = new GeneSurvival('fig');
	vis.setStatsDiv('#stats_div');
	
	%if 'recurrence' in locals():
		vis.x_label = 'Months since treatment started';
		vis.y_label = 'Recurrence free (%)';
	%end
	
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/groups_new.json';

	var test_group_selector = new GroupSelector('#test_groups', group_file,
		function() { downloaded(); });
	test_group_selector.title = 'Sample group';
	test_group_selector.default_group = 'Primary (untreated)';
	
	var downloaded = function() {
		if (!('sample_ids' in data && 'expr' in data && 
			'survival' in data && 'censored' in data &&
			test_group_selector.selected_samples.length > 0))
			return;
			
		// Pick test group samples.
		var whitelist = {};
		for (var s in test_group_selector.selected_samples)
			whitelist[test_group_selector.selected_samples[s]] = 1;
		
		var valid = function(val, idx) {
			return data.sample_ids[idx] in whitelist && 
				data.survival[idx] >= 0;
		};
		
		test_filt.expr = [];
		for (var e = 0; e < data.expr.length; e++)
			test_filt.expr[e] = data.expr[e].filter(valid);
		test_filt.survival = data.survival.filter(valid);
		test_filt.censored = data.censored.filter(valid);
		
		// Update the sample count display.
		$('#sample_count').html(
			test_filt.expr[0].length + ' / ' + data.expr[0].length);
		$('#censored_count').html(test_filt.censored.reduce(
			function(a, b) { return b == 1 ? a + 1 : a; }) + ' / ' +
			test_filt.censored.length);

		// Calculate differential splicing between the two groups.
		var e = data.selected_exon;
		
		var pad = 100;
		var ref_median = pv.median(data.expr[e], function(d) {
			if (d == 0) { pad = -pad; return pad; } else { return d; }
		});
		
		test_filt.diff = [];
		
		var exon_expr = test_filt.expr[e];
		for (var s = 0; s < exon_expr.length; s++) {
			if (exon_expr[s] == 0) continue;
			test_filt.diff.push(exon_expr[s] - ref_median);
		}
		
		vis.setData(test_filt.diff, test_filt.survival, test_filt.censored);
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
					downloaded();
					update_exon_list();
					return false;
				}
			})(e));
		}
	}

	var download_gene = function(gene_name) {
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
			downloaded();
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
		data.sample_ids = d['data']; downloaded();
	});
	
	$.getJSON(data_root + '/exon_splice/features.json', function(d) {
		$('#GeneSelect').omnicomplete(d['features'], function(val) {
			download_gene(val);
		});
	});
});
</script>

