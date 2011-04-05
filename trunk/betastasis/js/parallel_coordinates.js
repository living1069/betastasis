function ParallelCoordinates(div) {
	var w = 700, h = 400;
	var self = this;
	
	this.vis = new pv.Panel()
		.top(40).left(60)
		.width(w).height(h+100)
		.canvas(div);
	
	this.rule = this.vis.add(pv.Rule)
		.height(h);
		
	this.rule_label = this.rule.anchor('top').add(pv.Label)
		.top(-18).font('bold 10px sans-serif');
	this.rule_max = this.rule.anchor('top').add(pv.Label)
		.top(-6).font('bold 10px sans-serif');
	this.rule_min = this.rule.anchor('bottom').add(pv.Label)
		.top(h+5).font('bold 10px sans-serif');
	this.rule_corr = this.vis.add(pv.Label)
		.top(h+50).font('12px sans-serif').textAlign('center');
	
	// Create a panel for each sample.
	this.bg_line_panels = this.vis.add(pv.Panel)
		.height(h);
	this.bg_lines = this.bg_line_panels.add(pv.Line)
		.strokeStyle('#ccc')
		.lineWidth(1);
	
	this.line_panels = this.vis.add(pv.Panel)
		.height(h);
	this.lines = this.line_panels.add(pv.Line)
		.lineWidth(1);
	
	this.active_handle = this.vis.add(pv.Bar);
	this.handle = this.vis.add(pv.Bar);

	
	

	this.setData = function(dims, val) {
		var max = [], min = [], mean = [], stdev = [], pearson = [];
		for (var k = 0; k < val.length; k++) {
			max[k] = pv.max(val[k]);
			min[k] = pv.min(val[k]);
			var dx = (max[k] - min[k]) / 20;
			max[k] = Math.ceil((max[k] + dx) * 10) / 10;
			var new_min = Math.floor((min[k] - dx) * 10) / 10;
			if (min[k] < 0 || new_min >= 0)
				min[k] = new_min;
			
			mean[k] = pv.mean(val[k]);
			stdev[k] = pv.deviation(val[k]);
		}
		
		for (var d = 0; d < val.length - 1; d++) {
			pearson[d] = 0;
			for (var k = 0; k < val[d].length; k++) {
				pearson[d] += (val[d][k] - mean[d]) * (val[d+1][k] - mean[d+1]);
			}
			pearson[d] /= stdev[d] * stdev[d+1] * (val[d].length - 1);
			pearson[d] = Math.round(pearson[d] * 1000) / 1000;
		}
		
		var val = pv.transpose(val);
		
		this.dims = dims;
		this.val = val;
		
		var x = pv.Scale.linear(0, dims.length).range(0, w);
		var shade = pv.Scale.linear(min[0], max[0]).range(240, 0);
		
		var y = [], filter = [], active = [];
		for (k = 0; k < dims.length; k++) {
			y[k] = pv.Scale.linear(min[k], max[k]).range(0, h);
			filter[k] = { min: min[k], max: max[k] };
			active[k] = false;
		}
		
		this.rule.data(dims).left(function(d) { return x(this.index); });
		this.rule_label.text(function(d) { return dims[this.index]; });
		this.rule_max.text(function(d) { return '' + max[this.index]; });
		this.rule_min.text(function(d) { return '' + min[this.index]; });
		
		this.rule_corr.data(pearson)
			.left(function(d) { return x(this.index) + 70; })
			.text(function(d) { return sprintf('Pearson = %.3f', d); });
				
		this.bg_line_panels.data(val);
		this.bg_lines.data(dims)
			.left(function(t, d) { return x(this.index); })
			.bottom(function(t, d) { return y[this.index](d[this.index]);} );
		
		this.line_panels.data(val)
			.visible(function(d) {
				for (k = 0; k < dims.length; k++) {
					if (d[k] < filter[k].min || d[k] > filter[k].max)
						return false;
				}
				return true;
			});
		this.lines.data(dims)
			.strokeStyle(function(t, d) {
				return 'hsl(' + shade(d[0]) + ',50%,50%)';
			})
			.left(function(t, d) { return x(this.index); })
			.bottom(function(t, d) { return y[this.index](d[this.index]);} );
			
		function filterDim(d) {
			var t = this.index;
			filter[t].max = y[t].invert(h - d.y);
			filter[t].min = y[t].invert(h - d.y - d.dy);
			active[t] = true;
			self.active_handle.render();
			return false;
		}

		function unfilterDim(d) {
			var t = this.index;
			if (d.dy < 3) {
				filter[t].min = min[t];
				filter[t].max = max[t];
				//d.y = 0; d.dy = h;
				active[t] = false;
			}
			self.vis.render();
			return false;
		}
		
		this.active_handle.data(dims)
			.left(function(t) { return x(this.index) - 5; })
			.top(function(t) {
				return h - y[this.index](filter[this.index].max); })
			.width(10)
			.height(function(t) {
				return y[this.index](filter[this.index].max) -
					y[this.index](filter[this.index].min); })
			.strokeStyle('none')
			.fillStyle(function(t)
				{ return active[this.index] ? "rgba(255,0,0,.5)" : "none"; });

		this.handle
			.data(dims.map(function(dim) { return {y:0, dy:h, dim: dim}; }))
			.left(function(t) { return x(this.index) - 5; })
			.width(10).height(h)
			.strokeStyle('none')
			.fillStyle("hsla(0,0,80%,.3)")
			.cursor("crosshair")
			.event("mousedown", pv.Behavior.select())
			.event("select", filterDim)
			.event("selectend", unfilterDim);

		this.vis.render();
	}
	
	this.vis.render();
}

