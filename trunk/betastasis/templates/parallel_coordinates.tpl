%rebase base path=path, url=url, scripts='ui,parallel_coordinates,feature_selector'

%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Parallel coordinates</h3>
<div id="fig"></div>
</div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
Samples: <span id="sample_count"></span><br>
</p>
	
<div id="gene_select_div"></div><br>

<b>Data scale</b><br>
<div id="plot_scale">
	<input type="radio" id="log2_scale" name="plot_scale" /><label for="log2_scale">Log-2</label>
	<input type="radio" id="rank_scale" name="plot_scale" /><label for="rank_scale">Rank</label>
</div><br>



<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}<br>

</div>

<div id="sample_groups" class="sidebox"></div>
</div>

<script type="text/javascript">
var data_root = '{{data}}';
	
$(document).ready(function() {
	
	var group_file = '{{groups if 'groups' in locals() else ''}}';
	if (group_file == '')
		group_file = data_root + '/groups/groups_new.json';
	
	var vis = new ParallelCoordinates('fig');
	
	var group_selector = new GroupSelector('#sample_groups', group_file,
		function() { downloaded(); });
	
	var data = {};
	var dims = [];
	var values = [];
	
	var downloaded = function() {
		if (!('sample_ids' in data &&
			group_selector.selected_samples.length > 0))
			return;
		
		// Check whether we have downloaded everything.
		for (var k = 0; k < values.length; k++) {
			if (values[k].length == 0) return;
		}
		
		// Discard samples which do not have an expression value for one of
		// the features.
		var valid = pv.repeat([1], values[0].length);
		for (var d = 0; d < values.length; d++) {
			var ds = values[d];
			for (var s = 0; s < ds.length; s++) {
				if (ds[s] == 0) valid[s] = 0;
			}
		}
		
		var whitelist = {};
		for (var s in group_selector.selected_samples)
			whitelist[group_selector.selected_samples[s]] = 1;
		
		var filt = [];
		for (var d = 0; d < values.length; d++) {
			filt[d] = values[d].filter(function(val, idx) {
				return valid[idx] && data.sample_ids[idx] in whitelist;
			});
		}
		
		// Update the sample count display.
		$('#sample_count').html(filt[0].length + ' / ' + values[0].length);
		
		vis.setData(dims, filt);
	}
	
	var download_features = function(features) {
		values = [];
		
		for (var k = 0; k < features.length; k++) {
			dims[k] = '';
			values[k] = [];
			
			// The extra wrapper around the closure is there to capture
			// the variables "k" and "features[k]" for the closure.
			$.getJSON(data_root + '/expr/' + features[k][0].toLowerCase() + 
				'/' + features[k] + '.json', 
				(function(idx, feature) { return function(d) {
					dims[idx] = feature;
					values[idx] = d['data'];
					downloaded();
				};})(k, features[k]));
		}

	};
	
	var default_features = {{features if 'features' in locals() else '["TP53", "PTEN", "EGFR", "PDGFRA", "MYC"]'}};
	
	$('#plot_scale').buttonset();
	$('#log2_scale').attr('checked', true);
	$('#log2_scale').button('refresh');

	$("#log2_scale").button().click(function() {
		vis.scale = 'log-2'; downloaded(); });
	$("#rank_scale").button().click(function() {
		vis.scale = 'rank'; downloaded(); });
		
	$.getJSON(data_root + '/expr/features.json', function(d) {
		var feature_select = new FeatureSelector('#gene_select_div',
			default_features, d['features']);
		feature_select.callback = download_features;
	});
	
	$.getJSON(data_root + '/clinical/sample_id.json', function(d) {
		data.sample_ids = d['data']; downloaded();
	});
	
	download_features(default_features);
});
</script>

