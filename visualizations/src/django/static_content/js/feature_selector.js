function FeatureSelector(div, defaults) {
	markup = '<p><b>Selected features</b><br>\n';
	for (var k = 0; k < defaults.length; k++) {
		markup += '<input type="text" id="GeneSelect_' + k +
			'" value="' + defaults[k] + '" /><br>\n';
	}
	markup += '</p>\n';
	
	$(div).html(markup);
	
	var features = defaults;
	var self = this;
	
	for (var k = 0; k < defaults.length; k++) {
		$('#GeneSelect_' + k).quickselect({
			maxVisibleItems: 10, minChars: 1, matchMethod: 'quicksilver',
			autoSelectFirst: false, selectSingleMatch: false, data: features});
		
		$('#GeneSelect_' + k).change(function() {
			gs = $(this);
			// For some reason the $(this).val() hasn't yet updated when we
			// get to this function if the user selects a gene with the mouse.
			// So we use a tiny timeout here.
			setTimeout(function() {
				feature = gs.val();
				for (var f = 0; f < features.length; f++) {
					if (features[f].toLowerCase() == feature.toLowerCase())
						break;
				}
				
				if (f == features.length) {
					gs.css('background-color', '#f66');
					return;
				}
			
				gs.css('background-color', 'white');
				features[k] = feature;
				self.callback(features);
			}, 10);
		});
	}
}

