function FeatureMatrix(div) {
	var w = 22, h = 12;
	var feature_labelsize = 100;
	var legend_width = 100;
	this.sample_labelsize = 70;
	
	var vis, col_panels, matrix, col_labels, row_labels;
	
	this.selected = [];
	this.dragged_rows = [];
	this.highlighted_feature = -1;
	this.last_highlighted_feature = -1;
	
	self = this;
	
	function selectRows(d) {
		var start_row = Math.floor(d.y / h);
		var end_row = Math.floor((d.y + d.dy) / h);
		self.dragged_rows = [start_row, end_row];
		self.drag_bar.render();
		return false;
	}

	function endSelectRows(d) {
		if (d.dy == 0) {
			var row = Math.floor(d.y / h);
			self.selected[row] = 1 - self.selected[row];
		} else if (self.dragged_rows.length == 2) {
			for (var k = self.dragged_rows[0]; k <= self.dragged_rows[1]; k++) {
				self.selected[k] = 1 - self.selected[k];
			}
		}
		self.dragged_rows = [];
		self.selection_bar.render();
		self.drag_bar.render();
		return false;
	}
	
	// Data must be given as an array of arrays, where the inner arrays
	// represent rows.
	this.setData = function(data, cols, rows) {
		this.cols = cols;
		this.rows = rows;
		this.data = data;
		
		this.selected = pv.repeat([0], this.rows.length);
				
		if (vis == undefined) {
			vis = new pv.Panel().canvas(div);
			col_panels = vis.add(pv.Panel);
			matrix = col_panels.add(pv.Panel);
			col_labels = vis.add(pv.Label);
			row_labels = vis.add(pv.Label);
			this.selection_bar = vis.add(pv.Panel);
			this.drag_bar = vis.add(pv.Panel);
			this.crosshair_panel = vis.add(pv.Panel);
			
			this.legend_title = vis.add(pv.Label);
			this.legend = vis.add(pv.Dot);
			
			this.feature_buttons = vis.add(pv.Panel);
			this.feature_delete = vis.add(pv.Image)
				.width(w-4).height(w-4).url('http://viherjora.cs.tut.fi/betadata/delete.png');
		}
		
		this.updateData();
	}
	
	this.updateData = function() {
		if (vis == undefined) return;
		
		var show_sample_labels = (h >= 8);
		var sample_labelsize = show_sample_labels ? self.sample_labelsize : 0;
		
		vis.width(sample_labelsize + this.cols.length * w + 20 + legend_width)
			.height(feature_labelsize + this.rows.length * h + 30);

		col_panels.data(this.data)
			.top(feature_labelsize)
			.left(function() { return sample_labelsize + this.index * w; })
			.width(w);
		
		matrix.data(this.rows)
			.top(function() { return this.index * h; })
			.height(h)
			.fillStyle(function(t,d) { return self.cellColor(d[this.index]); })
			.strokeStyle("white")
			.lineWidth(1)
			.antialias(false);
		
		if (h < 3)
			matrix.strokeStyle('none');

		col_labels.data(this.cols)
			.top(feature_labelsize - 5)
			.left(function() { return sample_labelsize + this.index*w + w/2; })
			.textAngle(-Math.PI / 2)
			.textBaseline("middle");

		if (show_sample_labels) {
			row_labels.data(this.rows)
				.visible(true)
				.top(function()
					{ return feature_labelsize + this.index*h + h/2; })
				.left(sample_labelsize - 5)
				.textAlign('right').textBaseline('middle');
		} else {
			row_labels.visible(false);
		}
		
		this.selection_bar.data(this.selected)
			.left(sample_labelsize + this.cols.length * w + 8)
			.top(function() { return feature_labelsize + this.index * h; })
			.height(h)
			.width(3)
			.fillStyle(function(d) { return d ? 'black':'none'; })
			.strokeStyle("none")
			.antialias(false);
		
		this.drag_bar.data(this.selected)
			.left(sample_labelsize + this.cols.length * w + 3)
			.top(function() { return feature_labelsize + this.index * h; })
			.height(h)
			.width(3)
			.fillStyle(function() {
				if (self.dragged_rows.length == 2 &&
					this.index >= self.dragged_rows[0] &&
					this.index <= self.dragged_rows[1])
					return '#ccc';
				return 'none';
			})
			.strokeStyle("none")
			.antialias(false);

		this.crosshair_panel.data([{ y: 0, dy: 0 }])
			.top(feature_labelsize).left(sample_labelsize)
			.width(this.cols.length * w).height(this.rows.length * h)
			.fillStyle('hsla(0,0,0,0.001)')
			.cursor('crosshair')
			.event("mousedown", pv.Behavior.select())
			.event("select", selectRows)
			.event("selectend", endSelectRows);
		
		this.feature_buttons.data(this.cols)
			.left(function() { return sample_labelsize + this.index * w; })
			.top(0).width(w).height(feature_labelsize)
			.fillStyle(function() {
				return (this.index == self.highlighted_feature) ?
					'hsla(0,0,0,0.1)' : 'hsla(0,0,0,0.001)';
			})
			.event('mouseover', function() { 
				self.highlighted_feature = this.index;
				self.feature_buttons.render();
				self.feature_delete.render();
			})
			.event('mouseout', function() {
				self.last_highlighted_feature = self.highlighted_feature;
				self.highlighted_feature = -1;
				setTimeout(function() {
					self.feature_buttons.render();
					self.feature_delete.render();
				}, 10);
			})
			.event('click', function() { self.sortByFeature(this.index); });
		
		this.feature_delete
			.visible(function() { 
				return self.highlighted_feature != -1 && self.cols.length > 1;
			})
			.top(2).left(function() {
				return sample_labelsize + self.highlighted_feature * w + 2;
			})
			.event('mouseover', function() { 
				self.highlighted_feature = self.last_highlighted_feature;
			})
			.event('click', function() {
				var tmp = self.highlighted_feature;
				self.highlighted_feature = -1;
				self.featureClicked(tmp);
			});
		
		this.legend_title.data(['Legend'])
			.font('14px sans-serif')
			.left(sample_labelsize + this.cols.length * w + 22)
			.top(function() { return feature_labelsize + 14 + 12*this.index; });

		this.legend.data(['+1.5 .. inf', '+0.5 .. +1.5', '-0.5 .. +0.5', '-0.5 .. -1.5', '-1.5 .. -inf', '', 'Non-mutated', 'Mutated'])
			.left(sample_labelsize + this.cols.length * w + 30)
			.top(function() { return feature_labelsize + 22 + 12 * this.index; })
			.fillStyle(function(d) {
				var cell_vals = [ 2, 1, 0, -1, -2, -100, 0, 100 ];
				return self.cellColor(cell_vals[this.index]); })
			.strokeStyle(null)
		  .anchor('right').add(pv.Label)
			.text(function(d) { return d; });
		
		vis.render();
	}
	
	this.setColumnWidth = function(width) {
		w = width;
		this.updateData();
	}
	
	this.setRowHeight = function(height) {
		h = height;
		this.updateData();
	}
	
	this.selectAllRows = function() {
		for (var k = 0; k < this.selected.length; k++)
			self.selected[k] = 1;
		self.selection_bar.render();
	}
	
	this.deselectAllRows = function() {
		for (var k = 0; k < this.selected.length; k++)
			self.selected[k] = 0;
		self.selection_bar.render();
	}
	
	this.selectedRows = function() {
		var ret = '';
		for (var k = 0; k < this.selected.length; k++) {
			if (this.selected[k] == 1) {
				if (ret.length > 0) ret = ret + ', ';
				ret = ret + this.rows[k];
			}
		}
		if (ret.length == 0) {
			ret = 'No rows selected.';
		}
		return ret;
	}
}

