function GeneSurvival(div) {
	var w = 350, h = 200, sw = 200, sh = 100;
	var colors = ['#a55', '#55a'];
	
	this.vis = new pv.Panel()
		.width(w+sw+200).height(h+sh)
		.canvas(div);

	this.survbox = this.vis.add(pv.Panel)
		.width(w).height(h)
		.left(35).top(5);

	this.histbox = this.vis.add(pv.Panel)
		.width(sw).height(sh)
		.top(5).left(425);

	this.vis.render();
	
	var kaplan_meier = function(surv, censored) {
		var at_risk = surv.length;
		var prod = 1;
		var km = [];
		for (var k = 0; k < surv.length; k++) {
			if (censored[k] == 0) {
				prod *= ((at_risk - 1) / at_risk);
			}
			km[k] = prod;
			at_risk -= 1;
		}
		
		return pv.range(-1, surv.length).map(function(k) {
			if (k == -1) return { x: 0, y: 100 };
			return { x: surv[k] / 30, y: km[k] * 100 };
		});
	}
	
	this.updateSurvival = function() {
		if (this.orig_expr == undefined || this.survival == undefined)
			return;
		
		// Sort the expression values according to the sort order of the
		// survival times.
		if (this.expr == undefined) {
			this.expr = [];
			for (var k = 0; k < this.survival_order.length; k++) {
				this.expr[k] = this.orig_expr[this.survival_order[k]];
			}
		}
			
		var low = [], high = [], low_censored = [], high_censored = [];
		for (var k = 0; k < this.expr.length; k++) {
			if (this.expr[k] < this.expr_threshold) {
				low.push(this.survival[k]);
				low_censored.push(this.censored[k]);
			} else {
				high.push(this.survival[k]);
				high_censored.push(this.censored[k]);
			}
		}
		
		this.low_line.data(kaplan_meier(low, low_censored));
		this.high_line.data(kaplan_meier(high, high_censored));
		
		this.survbox.render();
	}
	
	this.cursorMoved = function() {
		this.cursor.left(this.histbox.mouse().x);
		this.cursor.render();
		
		this.expr_threshold = this.expr_min + this.histbox.mouse().x /
			sw * (this.expr_max - this.expr_min);
		this.updateSurvival();
	}
	
	this.updateHistogram = function() {
		var num_bins = 50;
		var bins = pv.range(this.expr_min, this.expr_max, 
			(this.expr_max - this.expr_min) / (num_bins - 1));
		var dist = pv.repeat([0], bins.length + 1);
		
		var bin = 0;
		for (var k = 0; k < this.sorted_expr.length; k++) {
			while (bins[bin] < this.sorted_expr[k]) {
				bin += 1;
			}
			dist[bin] += 1;
		}
		
		var dist_max = pv.max(dist);
		
		var sx = pv.Scale.linear(this.expr_min, this.expr_max).range(0, sw);
		var sy = pv.Scale.linear(0, dist_max).range(0, sh);
		
		if (this.histbox_xticks == undefined) {
			this.histbox_xticks = this.histbox.add(pv.Rule);
			this.histbox_yticks = this.histbox.add(pv.Rule);
			
			this.histogram = this.histbox.add(pv.Bar)
				.left(function() { return this.index * sw / num_bins; })
				.width(sw / num_bins)
				.bottom(0)
				.lineWidth(1).fillStyle('#aaa')
				.antialias(false);
			
			this.histbox.add(pv.Label)
				.data(["Expression threshold"])
				.left(sw/2).bottom(-35)
				.font("10pt Arial").textAlign("center");

				
			this.cursor = this.histbox.add(pv.Rule)
				.left(sw / 2).top(-4).bottom(-4)
				.strokeStyle("red");
			
			var foo = this;
			this.histbox.events("all").event("mousemove", function() {
				foo.cursorMoved() });
		}
		
		this.histbox_xticks.data(sx.ticks())
			.visible(function(d) { return d > 0; })
			.left(sx)
			.strokeStyle("#eee")
		  .add(pv.Rule)
			.bottom(-5)
			.height(5)
			.strokeStyle("#000")
		  .anchor("bottom").add(pv.Label)
			.text(sx.tickFormat);
			
		this.histbox_yticks.data(sy.ticks(5))
			.bottom(sy)
			.strokeStyle(function(d) { return d ? "#eee" : "#000"; })
		  .anchor("left").add(pv.Label)
			.text(sy.tickFormat);
		
		var foo = this;
		this.histogram.data(dist)
			.height(function(d) { return d / dist_max * sh; })
			.events('all').event('mousemove', function() {
				foo.cursorMoved() });
		
		this.histbox.render();
	}
		
	this.setExpression = function(d) {
		this.orig_expr = d['expr'];
		this.sorted_expr = this.orig_expr.slice(0).sort(pv.naturalOrder);
		this.expr_min = pv.min(this.orig_expr);
		this.expr_max = pv.max(this.orig_expr);
		this.expr_threshold = (this.expr_min + this.expr_max) / 2;
		
		// Settings this.expr to undefined signals to updateSurvival() that
		// the expression values should be resorted according to 
		this.expr = undefined;
		
		this.updateHistogram();
		this.updateSurvival();
	}

	this.setSurvival = function(d) {
		this.survival = d['survival'];
		this.censored = d['censored'];
		this.survival_order = [];
		
		ab = [];
		for (var i = 0; i < this.survival.length; i++) {
			ab[i] = [this.survival[i], this.censored[i], i];
		}
		ab.sort(function(left, right) {
			return left[0] < right[0] ? -1 : 1;
		});
			
		for (var i = 0; i < this.survival.length; i++) {
			this.survival[i] = ab[i][0];
			this.censored[i] = ab[i][1];
			this.survival_order[i] = ab[i][2];
		}
		
		var x = pv.Scale.linear(0, pv.max(this.survival) / 30).range(0, w),
			y = pv.Scale.linear(0, 100).range(0, h);
		
		if (this.survbox_xticks == undefined) {
			this.survbox_xticks = this.survbox.add(pv.Rule)
				.data(x.ticks())
				.visible(function(d) { return d > 0; })
				.left(x)
				.strokeStyle("#eee")
			  .add(pv.Rule)
				.bottom(-5)
				.height(5)
				.strokeStyle("#000")
			  .anchor("bottom").add(pv.Label)
				.text(x.tickFormat);

			this.survbox.add(pv.Rule)
				.data(y.ticks(5))
				.bottom(y)
				.strokeStyle(function(d) { return d ? "#eee" : "#000"; })
			  .anchor("left").add(pv.Label)
				.text(y.tickFormat);
			
			this.survbox.add(pv.Label)
				.data(["Percent survival"])
				.left(-20).bottom(h/2)
				.font("10pt Arial").textAlign("center")
				.textAngle(-Math.PI/2);

			this.survbox.add(pv.Label)
				.data(["Months since initial diagnosis"])
				.left(w/2).bottom(-35)
				.font("10pt Arial").textAlign("center");
				
			this.low_line = this.survbox.add(pv.Line)
				.interpolate("step-after")
				.left(function(d) { return x(d.x); })
				.bottom(function(d) { return y(d.y); })
				.strokeStyle(colors[1])
				.lineWidth(2);
				
			this.high_line = this.survbox.add(pv.Line)
				.interpolate("step-after")
				.left(function(d) { return x(d.x); })
				.bottom(function(d) { return y(d.y); })
				.strokeStyle(colors[0])
				.lineWidth(2);
				
			this.survbox.add(pv.Dot)
				.data(['High expression', 'Low expression'])
				.right(100)
				.top(function() { return 10 + this.index * 12; })
				.fillStyle(function(d) { return colors[this.index]; })
				.strokeStyle(null)
			  .anchor("right").add(pv.Label)
				.text(function(d) { return d; });
		}
		
		this.updateSurvival();
	}
}

