<h3>Aberration matrix</h3>
<div id="fig"></div>
</div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count">218</span><br>

<div id="row_height_value">Row height: 12</div>
<div style="padding-left:10px" id="slider"></div>
<br></p>


<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br><br>

<b>Add features</b><br>
Copy number: <input style="float:right; width:140px; height:18px" type="text" id="cna_select" value="" /><br>
Mutation: <input style="float:right; width:140px; height:18px" type="text" id="mutation_select" value="" /><p></p>

<button id="default_features_button">Default features</button>
<button id="clear_features_button">Clear all</button>
<br><br>

<b>Selection</b><br>
<button id="select_all_button">Select all</button>
<button id="select_none_button">Select none</button>
<button id="sample_export_button">Save</button>
<br><br>

<b>Sorting</b><br>
<button id="clear_sort">Clear sorting</button>
</div>

<div id="sample_groups" class="sidebox"></div>

</div>

<script type="text/javascript">
var data_root = '{{data_root}}';

$(document).ready(function() {
	
	var data = new Object;
	data.values = [];
	data.features = [];
	data.sort_by = '';
	data.sort_ascending = false;
	
	var default_features = '{{default_features if 'default_features' in locals() else ''}}';
	if (default_features == '')
		default_features = 'AR,TP53,PTEN,MYC,MIIP,ERBB2,RB1,MDM2,NCOA2,AKT1';
	default_features = default_features.split(',');
	
	$('#slider').slider({ min: 1, max: 20, value: 12,
		slide: function(event, ui) {
			$('#row_height_value').html('Row height: ' + ui.value);
		},
		change: function(event, ui) {
			vis.setRowHeight(ui.value);
		}
	});
	
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/groups_new.json';
	
	var group_selector = new GroupSelector('#sample_groups', group_file,
		function() { downloaded(); });

	var vis = new FeatureMatrix('fig');
	vis.setColumnWidth(20);
	vis.sample_labelsize = {{sample_labelsize if 'sample_labelsize' in locals() else '100'}};
	vis.cellColor = function(d) {
		if (d == 0) {
			return 'hsl(0,0%,70%)';
		} else if (d == -100) {
			return 'hsl(0,0%,100%)';
		} else if (d == 100) {
			return 'hsl(120,50%,70%)';
		} else if (d > 1.5) {
			return 'hsl(0,100%,70%)';
		} else if (d > 0) {
			return 'hsl(0,50%,70%)';
		} else if (d < -1.5) {
			return 'hsl(240,100%,70%)';
		} else if (d < 0) {
			return 'hsl(240,50%,70%)';
		}
	};
	
	var downloaded = function() {
		if (!('sample_ids' in data && data.values.length > 0 && 
			group_selector.selected_samples.length > 0))
			return;
			
		// Check whether we have downloaded everything.
		for (var k = 0; k < data.values.length; k++) {
			if (data.values[k].length == 0) return;
		}
		
		var whitelist = {};
		for (var s in group_selector.selected_samples)
			whitelist[group_selector.selected_samples[s]] = 1;
		
		var valid = function(val, idx) {
			return data.sample_ids[idx] in whitelist;
		}
		
		var filt = [];
		for (var f = 0; f < data.values.length; f++) {
			filt[f] = data.values[f].filter(valid);
		}
		
		var filt_sample_ids = data.sample_ids.filter(valid);
		
		// Sort by a specific column if we have to.
		var sort_idx = -1;
		for (var f = 0; f < data.features.length; f++) {
			if (data.features[f] == data.sort_by) {
				sort_idx = f; break;
			}
		}
		
		if (sort_idx != -1) {
			filt.push(filt_sample_ids);
			
			pv.transpose(filt);
			if (data.sort_ascending) {
				filt.sort(function(left, right) {
					return left[sort_idx] < right[sort_idx] ? -1 : 1;
				});
			} else { 
				filt.sort(function(left, right) {
					return left[sort_idx] > right[sort_idx] ? -1 : 1;
				});
			}
			pv.transpose(filt);
			
			filt_sample_ids = filt.pop();
		}
		
		// Update the sample count display.
		$('#sample_count').html(
			filt[0].length + ' / ' + data.sample_ids.length);
		
		vis.setData(filt, data.features, filt_sample_ids);
		$('#sidebar').css('left', vis.total_width);
	}
	
	var add_feature = function(feature) {
		var F = data.values.length;
		data.values[F] = [];
		data.features[F] = feature;
		
		var pos = feature.indexOf(' (CNA)');
		if (pos != -1) {
			feature = feature.slice(0, pos);
			data.features[F] = feature + ' (C)';
			$.getJSON(data_root + '/cna/' + feature[0].toLowerCase() + 
				'/' + feature + '.json', function(d) {
				data.values[F] = d['data'];
				downloaded();
			});
		}
		
		var pos = feature.indexOf(' (MUT)');
		if (pos != -1) {
			feature = feature.slice(0, pos);
			data.features[F] = feature + ' (M)';
			$.getJSON(data_root + '/mutation/' + feature[0].toLowerCase() + 
				'/' + feature + '.json', function(d) {
				data.values[F] = d['data'].map(
					function(v) { return v == 1 ? 100 : v; });
				downloaded();
			});
		}
	}
	
	vis.featureClicked = function(feature) {
		// Remove a clicked feature.
		if (data.features.length <= 1) return;
		
		var keep = function(val, idx) { return idx != feature; };
		data.values = data.values.filter(keep);
		data.features = data.features.filter(keep);
		downloaded();
	};
	
	vis.sortByFeature = function(feature_idx) {
		if (data.sort_by == data.features[feature_idx]) {
			data.sort_ascending = !data.sort_ascending;
		} else {
			data.sort_by = data.features[feature_idx];
			data.sort_ascending = false;
		}
		downloaded();
	};

	$('#sample_export_button').button().click(function() {
		alert(vis.selectedRows());
	});
	
	$('#select_all_button').button().click(function() { vis.selectAllRows(); });
	$('#select_none_button').button().click(
		function() { vis.deselectAllRows(); });
	
	$('#default_features_button').button().click(function() {
		data.values = [];
		data.features = [];
		for (var k = 0; k < default_features.length; k++)
			add_feature(default_features[k]);
	});
	
	$('#clear_features_button').button().click(function() {
		data.values = [];
		data.features = [];
		vis.setData([], [], []);
		$('#sample_count').html('0 / ' + data.sample_ids.length);

	});

	$('#clear_sort').button().click(function() {
		data.sort_by = ''; downloaded();
	});
	
	$.getJSON(data_root + '/cna/features.json', function(d) {
		for (var k = 0; k < default_features.length; k++)
			add_feature(default_features[k]);

		$('#cna_select').omnicomplete(d['features'], function(val) {
			add_feature(val + ' (CNA)');
			setTimeout(function() { $('#cna_select').val(''); }, 10);
		});
	});
	
	$.getJSON(data_root + '/mutation/features.json', function(d) {
		$('#mutation_select').omnicomplete(d['features'], function(val) {
			add_feature(val + ' (MUT)');
			setTimeout(function() { $('#mutation_select').val(''); }, 10);
		});
	});

	$.getJSON(data_root + '/clinical/sample_id.json', function(d) {
		data.sample_ids = d['data']; downloaded();
	});
});
</script>

