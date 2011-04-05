function MetadataTable(div) {
	var w = 24, h = 13, labelsize = 100;
	
	var vis, col_panels, matrix, col_labels, row_labels;
	
	this.selected = [];
	this.dragged_rows = [];
	
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
		this.selected[2] = 1;
		this.selected[5] = 1;
		this.selected[6] = 1;
		this.selected[7] = 1;
				
		if (vis == undefined) {
			vis = new pv.Panel();
			col_panels = vis.add(pv.Panel);
			matrix = col_panels.add(pv.Panel);
			col_labels = vis.add(pv.Label);
			row_labels = vis.add(pv.Label);
			this.selection_bar = vis.add(pv.Panel);
			this.drag_bar = vis.add(pv.Panel);
			this.crosshair_panel = vis.add(pv.Panel);
		}
		
		pv.transpose(this.data);
		this.updateData();
	}
	
	this.cellColor = function(d) {
		return d ? '#e99':'#aaa';
	}
	
	this.updateData = function() {
		if (vis == undefined) return;
		
		vis.width(labelsize + this.cols.length * w + 20)
			.height(labelsize + this.rows.length * h)
			.canvas(div);

		col_panels.data(this.data)
			.top(labelsize)
			.left(function() { return labelsize + this.index * w; })
			.width(w);
		
		matrix.data(this.rows)
			.top(function() { return this.index * h; })
			.height(h)
			.fillStyle(function(t,d) { return self.cellColor(d[this.index]); })
			.strokeStyle("white")
			.lineWidth(1)
			.antialias(false);

		col_labels.data(this.cols)
			.top(labelsize - 5)
			.left(function() { return labelsize + this.index * w + w / 2; })
			.textAngle(-Math.PI / 2)
			.textBaseline("middle");

		row_labels.data(this.rows)
			.top(function() { return labelsize + this.index * h + h / 2; })
			.left(labelsize - 5)
			.textAlign('right')
			.textBaseline('middle');
		
		this.selection_bar.data(this.selected)
			.left(labelsize + this.cols.length * w + 8)
			.top(function() { return labelsize + this.index * h; })
			.height(h)
			.width(3)
			.fillStyle(function(d) { return d ? 'black':'none'; })
			.strokeStyle("none")
			.antialias(false);
		
		this.drag_bar.data(this.selected)
			.left(labelsize + this.cols.length * w + 3)
			.top(function() { return labelsize + this.index * h; })
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
			.top(labelsize).left(labelsize)
			.width(this.cols.length * w).height(this.rows.length * h)
			.fillStyle('hsla(0,0,0,0.01)')
			.cursor('crosshair')
			.event("mousedown", pv.Behavior.select())
			.event("select", selectRows)
			.event("selectend", endSelectRows);
		
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

