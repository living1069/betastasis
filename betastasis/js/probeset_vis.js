function ProbesetVis(div) {
	var w = 600, h = 400;
	
	var exon_h = 10;
	var cds_h = 4;
	
	var vis, boxes, x_ticks, x_axis, tx_labels, tx_label_y = [];
	var zoom_panel;
	var selected_elem = -1, longest_tx = 0;
	
	var is_chrome = navigator.userAgent.toLowerCase().indexOf('chrome') > -1;
	var zoom_speed = is_chrome ? 1 : 1/12;
	
	page = self;
	self = this;
	
	/*
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
	*/
	
	this.setGene = function(data) {
		this.transcripts = [];
		this.tx_names = [];
		tx_label_y = [];
		for (var tx in data) {
			this.tx_names.push(tx);
			this.transcripts.push(data[tx]);
		}
		
		longest_tx = 0;
		for (var k = 0; k < this.transcripts.length; k++) {
			longest_tx = Math.max(this.transcripts[k].length,
				longest_tx);
		}
		
		x = pv.Scale.linear(0, longest_tx).range(0, w);
		
		var y_base = 0;

		this.vis_elems = [];
		for (var k = 0; k < this.tx_names.length; k++) {
			var tx = this.transcripts[k];
			
			var overlap = pv.repeat([0], tx.length);
			var highest_level = 0;
			
			for (var p = 0; p < tx.probe_pos.length; p++) {
				var probe = tx.probe_pos[p];
				var max_overlap = 0;
				for (var s = probe[0]; s <= probe[1]; s++) {
					if (overlap[s] > max_overlap)
						max_overlap = overlap[s];
				}
				
				var level = max_overlap + 1;
				if (level > highest_level)
					highest_level = level;
				
				for (var s = probe[0]; s <= probe[1]; s++) {
					overlap[s] = level;
				}
			}
			
			overlap = pv.repeat([0], tx.length);
			y_base += highest_level * 5;

			for (var p = 0; p < tx.probe_pos.length; p++) {
				var probe = tx.probe_pos[p];
				var max_overlap = 0;
				for (var s = probe[0]; s <= probe[1]; s++) {
					if (overlap[s] > max_overlap)
						max_overlap = overlap[s];
				}
				
				var level = max_overlap + 1;
				if (level > highest_level)
					highest_level = level;
				
				for (var s = probe[0]; s <= probe[1]; s++) {
					overlap[s] = level;
				}
				
				var color = 0.0;
				if (tx.probe_accepted[p] == 0)
					color = 0.7;
				
				this.vis_elems.push({ x: probe[0], y: y_base - 5 * level,
					w: probe[1] - probe[0], h: 2, color: color,
					tooltip: 'Sequence: ' + tx.probe_seq[p].toUpperCase() +
						'   \n' + 'Position: ' + probe[0] + ' .. ' + probe[1]});
			}
			
			for (var e = 0; e < tx.exons.length-1; e++) {
				var exon_start = tx.exon_pos[e];
				var exon_len = tx.exon_pos[e+1] - tx.exon_pos[e];
				this.vis_elems.push({ x: exon_start, y: y_base,
					w: exon_len-1, h: exon_h, color: 0.7,
					tooltip: this.tx_names[k] + ' exon ' + tx.exons[e] + 
						'   \nPosition: ' + exon_start + ' .. ' +
						(exon_start + exon_len - 1)});
			}
			
			this.vis_elems.push({ x: tx.exon_pos[tx.exons.length-1], y: y_base,
				w: tx.length - tx.exon_pos[tx.exons.length-1], h: exon_h,
				color: 0.7,
				tooltip: this.tx_names[k] + ' exon ' + 
					tx.exons[tx.exons.length-1] + 
					'   \nPosition: ' + tx.exon_pos[tx.exons.length-1] + 
					' .. ' + tx.length});
			
			if (tx.cds != undefined) {
				this.vis_elems.push({ x: tx.cds[0], y: y_base + 3,
					w: tx.cds[1] - tx.cds[0], h: exon_h - 6, color: 0.3,
					tooltip: this.tx_names[k] + ' CDS' +
						'   \nPosition: ' + tx.cds[0] + ' .. ' + tx.cds[1] });
			}

			
			tx_label_y[k] = y_base + exon_h + 2;
			y_base += exon_h + 20;
		}
		
		if (vis == undefined) {
			vis = new pv.Panel();
			zoom_panel = vis.add(pv.Panel);
			ivis = vis.add(pv.Panel);
			boxes = ivis.add(pv.Panel);
			tx_labels = ivis.add(pv.Label);
			tx_label_links = ivis.add(pv.Panel);
			x_axis = ivis.add(pv.Rule);
			x_ticks = ivis.add(pv.Rule);
			
			/*
			matrix = col_panels.add(pv.Panel);
						row_labels = vis.add(pv.Label);
			this.selection_bar = vis.add(pv.Panel);
			this.drag_bar = vis.add(pv.Panel);
			*/
		}
		
		this.updateData();
	}
	
	this.updateData = function() {
		if (vis == undefined) return;
		
		vis.width(w + 3)
			.height(tx_label_y[tx_label_y.length-1] + 40)
			.canvas(div);
		
		ivis.width(w).left(3);
		
		zoom_panel.cursor('crosshair')
			.events("all")
			.event("mousedown", pv.Behavior.pan().bound(true))
			.event("mousewheel", pv.Behavior.zoom(zoom_speed).bound(true))
			.event("pan", transform)
			.event("zoom", transform);
		
		boxes.data(this.vis_elems)
			.top(function(d) { return d.y; })
			.left(function(d) { return x(d.x); })
			.width(function(d) { return x(d.x + d.w) - x(d.x); })
			.height(function(d) { return d.h; })
			.fillStyle(function(d) {
				if (this.index == self.selected_elem) {
					return '#66f';
				} else {
					return 'hsl(0,0,' + Math.round(d.color * 100) + '%)';
				}
			})
			.title(function(d) { return d.tooltip; })
			.cursor('pointer')
			.event('mouseover', function() {
				self.selected_elem = this.index; boxes.render(); })
			.event('mouseout', function() {
				self.selected_elem = -1; boxes.render(); });

		tx_labels.data(this.tx_names)
			.top(function() { return tx_label_y[this.index]; })
			.left(-3)
			.textAlign('left')
			.textBaseline('top');
		
		tx_label_links.data(this.tx_names)
			.top(function() { return tx_label_y[this.index] - 2; })
			.left(-5).width(90).height(14)
			.fillStyle('rgba(1.0,0,0,0.001)')
			.cursor('pointer')
			.event('click', function(d) {
				window.open('http://www.ncbi.nlm.nih.gov/nuccore?term=' + d, '_blank');
			});

		x_axis.data([0]).bottom(20);
		
		x_ticks
			.data(function() {
				var x_min = pv.min(x.domain()), x_max = pv.max(x.domain());
				if (x_max - x_min < w / 20) {
					return pv.range(Math.round(x_min), Math.round(x_max));
				} else {
					return x.ticks();
				}
			})
			.visible(true)
			.left(x).bottom(15).height(5)
		  .anchor('bottom').add(pv.Label)
			.text(function(d) { return d; });
		
		function transform() {
			var t = this.transform().invert();
			x.domain(t.x / w * longest_tx, (t.k + t.x / w) * longest_tx);
			ivis.render();
		}
		
		/*
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
		*/
		
		vis.render();
	}
	
	/*
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
	*/
}

