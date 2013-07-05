%rebase base path=path, url=url, scripts='ui,feature_matrix'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Aberration matrix</h3>
<div id="fig"></div>
</div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count"></span><br>

<div id="row_height_value">Row height: 12</div>
<div style="padding-left:10px" id="slider"></div>
<div>Expression Treshold</div>
<div style="padding-left:10px" id="expr_slider"></div>
<div>CNA Treshold</div>
<div style="padding-left:10px" id="cna_slider"></div>
<br></p>

<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br><br>

<b>Add features</b><br>
Expression: <input style="float:right; widhth:140px; height:18px" type="text" id="expr_select" value="" /><br>
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
<br><br>

<b>Export</b><br>
<button id="export_table">Table</button>
<button id="export_svg">SVG</button>
<br>

</div>

<div id="sample_groups" class="sidebox"></div>

</div>

<script type="text/javascript">
var data_root = '{{data}}';
var treshhold = {};
treshhold['cna'] = [0,-100,100,1.5,0.5,-0.5,-1.5];
treshhold['expr'] = [0,-100,100, 1.5, 0.5, -0.5, -1.5];
treshhold['mut'] = [0,1,2,3,4];
treshhold['legend'] = ['+1.5 .. inf', '+0.5 .. +1.5', '-0.5 .. +0.5', '-0.5 .. -1.5', '-1.5 .. -inf', '', '', '+1.5 .. inf', '+0.5 .. +1.5', '-0.5 .. +0.5', '-0.5 .. -1.5', '-1.5 .. -inf', '', '', 'Non-mutated', 'Heterozygote non-synonymous mutation', 'Homozygote non-synonymous mutation', 'Heterozygote silent mutation', 'Homozygote silent mutation']
var vis;



$(document).ready(function() {
	
	var data = new Object;
	var filt = {};
	data.values = [];
	data.features = [];
	data.sort_by = '';
	data.sort_ascending = false;
	
	var default_features = '{{default_features if 'default_features' in locals() else ''}}';
	if (default_features == '')
		default_features = 'AR (CNA),TP53 (CNA),PTEN (CNA),MYC (CNA),MIIP (CNA),ERBB2 (CNA),RB1 (CNA),MDM2 (CNA),NCOA2 (CNA),AKT1 (CNA)';
	default_features = default_features.split(',');
	
	$('#slider').slider({ min: 1, max: 20, value: 12,
		slide: function(event, ui) {
			$('#row_height_value').html('Row height: ' + ui.value);
		},
		change: function(event, ui) {
			vis.setRowHeight(ui.value);
		}
	});

	initSlider('expr');
	initSlider('cna');
	function initSlider (id) {

		var index = (id === 'cna') ? 7 : 0;
		var init_val = (id === 'cna') ? [ -150, -50, 50, 150 ] : [ -150, -50, 50, 150];
	$( '#' + id+'_slider' ).slider({
      range: false,
      min: -400,
      max: 400,
      values: init_val,
      change: function( event, ui ) {
        
        treshhold[id][3] = ui.values[3]/100;
        treshhold[id][4] = ui.values[2]/100;
        treshhold[id][5] = ui.values[1]/100;
        treshhold[id][6] = ui.values[0]/100;
        treshhold['legend'][index] = ((ui.values[3]/100.0 > 0.0) ? '+':'') + (ui.values[3]/100.0).toFixed(1) + ' .. inf';
        treshhold['legend'][index + 1] = ((ui.values[2]/100.0 > 0) ? '+':'') + (ui.values[2]/100.0).toFixed(1) + ' .. ' + ((ui.values[3]/100.0 > 0) ? '+':'') + (ui.values[3]/100.0).toFixed(1);
        treshhold['legend'][index + 2] = ((ui.values[1]/100 > 0) ? '+':'') + (ui.values[1]/100).toFixed(1) + ' .. ' + ((ui.values[2]/100 > 0) ? '+':'') + (ui.values[2]/100).toFixed(1);
        treshhold['legend'][index + 3] = ((ui.values[1]/100 > 0) ? '+':'') + (ui.values[1]/100).toFixed(1) + ' .. ' + ((ui.values[0]/100 > 0) ? '+':'') + (ui.values[0]/100).toFixed(1);
        treshhold['legend'][index + 4] = ((ui.values[0]/100 > 0) ? '+':'') + (ui.values[0]/100).toFixed(1) + ' .. -inf';
        vis.updateData();
      }
    });
  	handle = $('#' + id + '_slider' +' A.ui-slider-handle');        
    handle.eq(0).addClass('first-handle');        
    handle.eq(1).addClass('second-handle');
    handle.eq(2).addClass('third-handle');
    handle.eq(3).addClass('fourth-handle');
	}



	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/groups_new.json';
	
	var group_selector = new GroupSelector('#sample_groups', group_file,
		function() { downloaded(); });

	 vis = new FeatureMatrix('fig');
	vis.setColumnWidth(20);
	vis.sample_labelsize = {{sample_labelsize if 'sample_labelsize' in locals() else '100'}};
	vis.cellColor = function(d, id) {
		
		if(id !== undefined){
			if(id.search('(C)') !== -1){
				if (d >= treshhold['cna'][3]) {
					return 'hsl(0,100%,70%)';
				} else if (d >= treshhold['cna'][4]) {
					return 'hsl(0,50%,70%)';
				} else if (d <= treshhold['cna'][6]) {
					return 'hsl(240,100%,70%)';
				} else if (d <= treshhold['cna'][5]) {
					return 'hsl(240,50%,70%)';
				} else if(d <= treshhold['cna'][4] && d >= treshhold['cna'][5]){
					return 'hsl(0,0%,70%)';
				}
			}
			if(id.search('(E)') !== -1){
				if (d >= treshhold['expr'][3]) {
					return 'hsl(0,100%,70%)';
				} else if (d >= treshhold['expr'][4]) {
					return 'hsl(0,50%,70%)';
				} else if (d <= treshhold['expr'][6]) {
					return 'hsl(240,100%,70%)';
				} else if (d <= treshhold['expr'][5]) {
					return 'hsl(240,50%,70%)';
				} else if(d <= treshhold['expr'][4] && d >= treshhold['expr'][5]){
					return 'hsl(0,0%,70%)';
				}
			}
			if(id.search('(M)') !== -1){
				var color_str ="";
				switch(d)
				{
					case 0:
						color_str = 'hsl(0,0%,70%)';
						break;
					case 1:
						color_str = 'hsl(120,100%, 70%)';
						break;
					case 2:
						color_str = 'hsl(120,100%, 40%)';
						break;
					case 3:
						color_str = 'hsl(60,100%, 70%)';
						break;
					case 4:
						color_str = 'hsl(60,100%, 40%)';
						break;
					default:
						color_str = 'hsl(0,0%, 100%)';
				}
				return color_str;
			}
		}
		else{ 
			var color_str ="";
				switch(d)
				{
					case 0:
						color_str = 'hsl(0,0%,70%)';
						break;
					case 1:
						color_str = 'hsl(0,50%,70%)';
						break;
					case 2:
						color_str = 'hsl(0,100%,70%)';
						break;
					case -1:
						color_str = 'hsl(240,50%,70%)';
						break;
					case -2:
						color_str = 'hsl(240,100%,70%)';
						break;
					case 10:
						color_str = 'hsl(120,100%, 70%)';
						break;
					case 20:
						color_str = 'hsl(120,100%, 40%)';
						break;
					case 30:
						color_str = 'hsl(60,100%, 70%)';
						break;
					case 40:
						color_str = 'hsl(60,100%, 40%)';
						break;
					default:
						color_str = 'hsla(0,0%, 0%, 0)';
				}
				return color_str;
				
		}
		return 'hsl(0,0%,0%)';
	};
	
	var export_table = function() {
		var tabular = 'Sample';
		for (var k = 0; k < data.features.length; k++) {
			tabular += '\t' + data.features[k];
		}
		tabular += '\n';
		
		for (var s = 0; s < filt.rows.length; s++) {
			tabular += filt.rows[s];
			for (var k = 0; k < data.features.length; k++) {
				var pos = data.features[k].indexOf(' (C)');
				if (pos != -1) {
					var d = filt.data[k][s];
					if (d == 0) { tabular += '\tNeutral'; }
					else if (d > 0) { tabular += '\tAmplified'; }
					else if (d < 0) { tabular += '\tDeleted'; }
					else { tabular += '\t'; }
				} else {
					tabular += '\t';
				}
			}
			tabular += '\n';
		}
		
		download_data_uri(tabular, 'feature_matrix.txt', 
			'data:application/octet-stream');
	}
	
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
		
		filt.data = [];
		for (var f = 0; f < data.values.length; f++) {
			filt.data[f] = data.values[f].filter(valid);
		}
		
		filt.rows = data.sample_ids.filter(valid);
		
		// Sort by a specific column if we have to.
		var sort_idx = -1;
		for (var f = 0; f < data.features.length; f++) {
			if (data.features[f] == data.sort_by) {
				sort_idx = f; break;
			}
		}
		
		if (sort_idx != -1) {
			filt.data.push(filt.rows);
			
			pv.transpose(filt.data);
			if (data.sort_ascending) {
				filt.data.sort(function(left, right) {
					return left[sort_idx] < right[sort_idx] ? -1 : 1;
				});
			} else { 
				filt.data.sort(function(left, right) {
					return left[sort_idx] > right[sort_idx] ? -1 : 1;
				});
			}
			pv.transpose(filt.data);
			
			filt.rows = filt.data.pop();
		}
		
		// Update the sample count display.
		$('#sample_count').html(
			filt.data[0].length + ' / ' + data.sample_ids.length);
		
		vis.setData(filt, data.features);
		$('#sidebar').css('left', vis.total_width);
	}
	
	var add_feature = function(feature) {
		var F = data.values.length;
		data.values[F] = [];
		data.features[F] = feature;
		
		var pos = feature.indexOf(' (EXP)');
		if (pos != -1) {
			feature = feature.slice(0, pos);
			data.features[F] = feature + ' (E)';
			$.getJSON(data_root + '/diff_expr/' + feature[0].toLowerCase() + 
				'/' + feature + '.json', function(d) {
				data.values[F] = d['data'].map(function(x) {
					if (x == -100) return -100;
					else if (x > 3) return 2;
					else if (x > 1) return 1;
					else if (x < -3) return -2;
					else if (x < -1) return -1;
					else return 0;
				});
				downloaded();
			});
		}
		
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
	
	$('#export_table').button().click(export_table);
	$('#export_svg').button().click(function()
		{ export_svg('feature_matrix.svg'); });
	
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
		vis.setData({data: [], rows: []}, []);
		$('#sample_count').html('0 / ' + data.sample_ids.length);

	});

	$('#clear_sort').button().click(function() {
		data.sort_by = ''; downloaded();
	});
	
	$.getJSON(data_root + '/diff_expr/features.json', function(d) {
		$('#expr_select').omnicomplete(d['features'], function(val) {
			add_feature(val + ' (EXP)');
			setTimeout(function() { $('#expr_select').val(''); }, 10);
		});
		
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

