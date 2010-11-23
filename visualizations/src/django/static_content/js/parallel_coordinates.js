function ParallelCoordinates(div) {
	var w = 700, h = 400;
	
	var vis = new pv.Panel()
		.top(40).left(40)
		.width(w).height(h+20)
		.canvas(div);

	this.vis = vis;
	
	this.setData = function(dims, val) {
		max = []; min = [];
		for (k = 0; k < val.length; k++) {
			max[k] = pv.max(val[k]);
			min[k] = pv.min(val[k]);
			var dx = (max[k] - min[k]) / 20;
			max[k] = Math.ceil((max[k] + dx) * 10) / 10;
			var new_min = Math.floor((min[k] - dx) * 10) / 10;
			if (min[k] < 0 || new_min >= 0)
				min[k] = new_min;
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
		
		var rule = this.vis.add(pv.Rule)
			.data(dims).height(h)
			.left(function(d) { return x(this.index); });
	
		rule.anchor("top").add(pv.Label)
			.top(-18).font("bold 10px sans-serif")
			.text(function() { return dims[this.index]; });
		rule.anchor("top").add(pv.Label)
			.top(-6).font("bold 10px sans-serif")
			.text(function() { return '' + max[this.index]; });
		rule.anchor("bottom").add(pv.Label)
			.top(h+5).font("bold 10px sans-serif")
			.text(function() { return '' + min[this.index]; });
		
		this.bg_lines = this.vis.add(pv.Panel)
			.data(val)  // Create a panel for each sample
			.height(h)
		  .add(pv.Line)
			.data(dims)
			.left(function(t, d) { return x(this.index); })
			.bottom(function(t, d) { return y[this.index](d[this.index]);} )
			.strokeStyle('#ccc')
			.lineWidth(1);
		
		this.lines = this.vis.add(pv.Panel)
			.data(val)  // Create a panel for each sample
			.height(h)
			.visible(function(d) {
				for (k = 0; k < dims.length; k++) {
					if (d[k] < filter[k].min || d[k] > filter[k].max)
						return false;
				}
				return true;
			})

		  .add(pv.Line)
			.data(dims)
			.left(function(t, d) { return x(this.index); })
			.bottom(function(t, d) { return y[this.index](d[this.index]);} )
			.strokeStyle(function(t, d) {
				return 'hsl(' + shade(d[0]) + ',50%,50%)';
			})
			.lineWidth(1);
			
		function filterDim(d) {
			var t = this.index;
			filter[t].max = y[t].invert(h - d.y);
			filter[t].min = y[t].invert(h - d.y - d.dy);
			active[t] = true;
			active_handle.render();
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
			vis.render();
			return false;
		}
		
		var active_handle = this.vis.add(pv.Bar)
			.data(dims)
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

		var handle = this.vis.add(pv.Bar)
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

