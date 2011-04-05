jQuery.fn.omnicomplete = function(data, select) {
	var textbox = $(this[0]);
	var ci_data = data.slice(0);
	for (var k = 0; k < ci_data.length; k++) {
		ci_data[k] = ci_data[k].toUpperCase();
	}
	
	textbox.autocomplete({
		minLength: 0, delay: 0,
		source: function(req, add) {
			var suggestions = [];
			var best = [];
			
			req = req.term.toUpperCase();
			if (req.length < 1) { add([]); return; }
			
			for (var i = 0; i < data.length; i++) {
				var pos = ci_data[i].indexOf(req);
				if (pos != -1) {
					if (pos == 0) {
						best.push(data[i]);
						if (best.length == 20) break;
					} else {
						suggestions.push(data[i]);
					}
				}
			}
			suggestions = best.concat(suggestions);
			add(suggestions.slice(0, 20));
		},
		open: function() {
			var menu = $(this).data("autocomplete").menu;
			menu.activate( $.Event({ type: "mouseenter" }), 
				menu.element.children().first() );
		},
		search: function(event, ui) {
			if (ci_data.indexOf($(this).val().toUpperCase()) == -1) {
				textbox.css('background-color', '#faa');
			} else {
				textbox.css('background-color', 'white');
			}
		},
		select: function(event, ui) {
			textbox.css('background-color', 'white');
			select(ui.item.value);
		},
		change: function(event, ui) {
			select($(this).val());
		}
	});
}

