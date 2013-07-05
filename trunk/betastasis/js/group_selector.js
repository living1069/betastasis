function GroupSelector(div, json_file, select) {
	
	var self = this;
	
	self.title = 'Sample groups';
	self.ready = false;
	self.groups = [];
	self.selected = '';
	self.selected_samples = [];
	self.default_group = undefined;
	
	var render = function() {
		markup = '<p><b>' + self.title + '</b><br><ul>\n';
		var k = 0;
		for (var g in self.groups) {
			var group_name = g;
			if (g[0] == '>') {
				markup += '<li class="nested">';
				group_name = g.slice(1);
			} else {
				markup += '<li>';
			}
			
			if (group_name[0] == '#') {
				markup += '<a>' + group_name.slice(1) + '</a></li>';
			} else {
				markup += '<a id="' + div.slice(1) + '_sample_group_' + k + '" href="">';
				if (g == self.selected) markup += '<b>';
				markup += group_name;
				if (g == self.selected) markup += '</b>';
				markup += '</a></li>\n';
			}
			
			k++;
		}
		markup += '</ul></p>\n';
		$(div).html(markup);
		
		var k = 0;
		for (var g in self.groups) {
			if (g[0] != '#' && g[1] != '#') {
				$('#' + div.slice(1) + '_sample_group_' + k).click((function(group_name) {
					return function() {
						if (self.selected == group_name) return false;
						self.selected = group_name;
						self.selected_samples = self.groups[group_name];
						try {
							select();
						} catch(err) {
							alert(err.name + '\n' + err.message + '\n');
						}
						render();
						return false;
					}
				})(g));
			}
			k++;
		}
	}
	
	render();
	
	$.getJSON(json_file, function(d) {
		self.groups = d;
		for (var g in self.groups) { self.selected = g; break; }
		if (self.default_group != undefined) {
			for (var g in self.groups) {
				if (g == self.default_group) { self.selected = g; break; }
			}
		}
		self.selected_samples = self.groups[self.selected];
		
		self.ready = true;
		
		select();
		render();
	});
}

