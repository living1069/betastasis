function FeatureSelector(div, defaults, features) {
	markup = '<p><b>Selected features</b><br>\n';
	for (var k = 0; k < defaults.length; k++) {
		markup += '<input type="text" id="GeneSelect_' + k +
			'" value="' + defaults[k] + '" /><br>\n';
	}
	markup += '</p>\n';
	
	$(div).html(markup);
	
	var self = this;
	this.features = defaults;
	
	for (var k = 0; k < defaults.length; k++) {
		$('#GeneSelect_' + k).omnicomplete(features, (function(idx) {
			return function(val) {
				self.features[idx] = val;
				self.callback(self.features);
			};
		})(k));
	}
}

