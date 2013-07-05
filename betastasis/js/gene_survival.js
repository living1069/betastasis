function GeneSurvival(div) {
	var w = 350, h = 200, sw = 200, sh = 100;
	var colors = ['#a55', '#55a'];
	var self = this;
	
	self.x_label = 'Months since initial diagnosis';
	self.y_label = 'Percent survival';
	
	this.vis = new pv.Panel()
		.width(w+sw+200).height(h+50)
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
		
		var ret = [{ x: 0, y: 100 }];
		
		for (var k = 0; k < surv.length; k++) {
			if (censored[k] == 0) {
				prod *= ((at_risk - 1) / at_risk);
			}
			at_risk -= 1;
			
			ret[k+1] = { x: surv[k] / 30, y: prod * 100,
				censored: censored[k] };
		}
		return ret;
	}
	
	this.updateSurvival = function() {
		if (this.orig_expr == undefined || this.survival == undefined)
			return;
		
		var x = pv.Scale.linear(0, pv.max(this.survival) / 30).range(0, w);
		
		this.survbox_xticks.data(x.ticks())
			.visible(function(d) { return d > 0; })
			.left(x);
		
		this.survbox_xticklabels
			.text(x.tickFormat);
		
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
		
		this.low_censored.data(kaplan_meier(low, low_censored));
		this.high_censored.data(kaplan_meier(high, high_censored));
		
		this.survbox.render();
		
		if (this.stats_div != undefined) {
			$(this.stats_div).html(
				'<button id="stats_button">Calculate</button><br>');
			
			$('#stats_button').button();
			$('#stats_button').click(function() {
				pval = self.calculateP();
				if (pval > 0.001) {
					pval = pval.toPrecision(3);
				} else if (pval != 0) {
					pval = pval.toExponential(2);
				}
				$(self.stats_div).html("P-value (logrank test): " + pval);
			});
		}
	}
	
	this.setExprThreshold = function(val) {
		this.cursor.left(sw * (val - this.expr_min) / 
			(this.expr_max - this.expr_min));
		this.cursor.render();
		this.expr_threshold = val;
		this.updateSurvival();
	}
	
	this.cursorMoved = function() {
		this.cursor.left(this.histbox.mouse().x);
		this.cursor.render();
		
		this.expr_threshold = this.expr_min + this.histbox.mouse().x /
			sw * (this.expr_max - this.expr_min);
		this.updateSurvival();
	}
	
	// www.nmr.mgh.harvard.edu/Neural_Systems_Group/strang/python/stats.py
	var normcdf = function(z) {
		var Z_MAX = 6.0;    // maximum meaningful z-value
		if (z == 0.0) {
			var x = 0.0;
		} else {
			var y = 0.5 * Math.abs(z);
			if (y >= 0.5 * Z_MAX) {
				var x = 1.0;
			} else if (y < 1.0) {
				var w = y*y;
				var x = ((((((((0.000124818987 * w
							-0.001075204047) * w +0.005198775019) * w
						  -0.019198292004) * w +0.059054035642) * w
						-0.151968751364) * w +0.319152932694) * w
					  -0.531923007300) * w +0.797884560593) * y * 2.0;
			} else {
				y = y - 2.0;
				var x = (((((((((((((-0.000045255659 * y
								 +0.000152529290) * y -0.000019538132) * y
							   -0.000676904986) * y +0.001390604284) * y
							 -0.000794620820) * y -0.002034254874) * y
						   +0.006549791214) * y -0.010557625006) * y
						 +0.011630447319) * y -0.009279453341) * y
					   +0.005353579108) * y -0.002141268741) * y
					 +0.000535310849) * y +0.999936657524;
			}
		}
		
		if (z > 0.0) {
			var prob = (x+1.0)*0.5;
		} else {
			var prob = (1.0-x)*0.5;
		}
		return prob;
	}
	
	this.calculateP = function() {
		if (this.expr == undefined || this.survival == undefined)
			return 0;
		
		surv = this.survival;
		censored = this.censored;
		
		var LOW = 0, HIGH = 1;
		
		var times = [];
		var idx = 0;
		var at_risk = [ [0], [0] ];
		var deaths = [ [0, 0], [0, 0] ];
		
		for (var k = 0; k < this.expr.length; k++) {
			if (this.expr[k] < this.expr_threshold) {
				at_risk[LOW][0] += 1;
			} else {
				at_risk[HIGH][0] += 1;
			}
		}

		times[0] = surv[k];
		at_risk[LOW][1] = at_risk[LOW][0];
		at_risk[HIGH][1] = at_risk[HIGH][0];
		
		for (var k = 0; k < surv.length; k++) {
			if (k > 0 && surv[k] != surv[k-1]) {
				idx += 1;
				times[idx] = surv[k];
				
				at_risk[LOW][idx+1] = at_risk[LOW][idx];
				at_risk[HIGH][idx+1] = at_risk[HIGH][idx];
				
				deaths[LOW][idx+1] = 0;
				deaths[HIGH][idx+1] = 0;
			}
			
			group = HIGH;
			if (this.expr[k] < this.expr_threshold) {
				group = LOW;
			}
			
			at_risk[group][idx+1] -= 1;
			if (!censored[k]) {
				deaths[group][idx] += 1;
			}
		}
		
		var a = 0.0, b = 0.0;
		for (var t = 0; t < times.length - 1; t++) {
			// See http://en.wikipedia.org/wiki/Log_rank_test
			Oj = deaths[LOW][t] + deaths[HIGH][t];
			Nj = at_risk[LOW][t] + at_risk[HIGH][t];
			a += deaths[LOW][t] - Oj * at_risk[LOW][t] / Nj;
			b += Oj * (at_risk[LOW][t] / Nj) * 
				(1 - at_risk[LOW][t] / Nj) * (Nj - Oj) / (Nj - 1);
			//alert([Oj, Nj, a, b]);
		}
		
		var z = a / Math.sqrt(b);
		if (z < 0) {
			return 2 * normcdf(z);
		} else {
			return 2 * (1 - normcdf(z));
		}
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
		
		var dist_max = Math.max(pv.max(dist), 5);
		
		var sx = pv.Scale.linear(this.expr_min, this.expr_max).range(0, sw);
		var sy = pv.Scale.linear(0, dist_max).range(0, sh);
		
		if (self.histbox_xticks == undefined) {
			self.histbox_xticks = self.histbox.add(pv.Rule);
			self.histbox_yticks = self.histbox.add(pv.Rule);
			
			self.histbox_xtick_labels = self.histbox_xticks.anchor("bottom")
				.add(pv.Label);
			self.histbox_ytick_labels = self.histbox_yticks.anchor("left")
				.add(pv.Label);
			
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
			
			this.histbox.events("all").event("mousemove", function() {
				self.cursorMoved() });
		}
		
		this.histbox_xticks.data(sx.ticks(5))
			.left(sx)
			.strokeStyle("#eee")
		  .add(pv.Rule)
			.bottom(-5)
			.height(5)
			.strokeStyle("#000");
			
		this.histbox_yticks.data(sy.ticks(5))
			.bottom(sy)
			.strokeStyle(function(d) { return d ? "#eee" : "#000"; });
		
		this.histbox_xtick_labels.text(sx.tickFormat);
		this.histbox_ytick_labels.text(sy.tickFormat);
		
		this.histogram.data(dist)
			.height(function(d) { return d / dist_max * sh; })
			.events('all').event('mousemove', function() {
				self.cursorMoved() });
		
		this.histbox.render();
	}
	
	this.setData = function(expr, survival, censored) {
		this.orig_expr = expr.slice(0);
		this.survival = survival.slice(0);
		this.censored = censored.slice(0);
		
		this.sorted_expr = this.orig_expr.slice(0).sort(pv.naturalOrder);
		this.expr_min = pv.min(this.orig_expr);
		this.expr_max = pv.max(this.orig_expr);
		this.expr_threshold = (this.expr_min + this.expr_max) / 2;
		
		// Settings this.expr to undefined signals to updateSurvival() that
		// the expression values should be resorted according to the same order
		// as the survival times.
		this.expr = undefined;
		
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
				.bottom(-5)
				.height(5)
				.strokeStyle("#000")
		
			this.survbox_xticklabels = this.survbox_xticks
			  .anchor("bottom").add(pv.Label);

			this.survbox.add(pv.Rule)
				.data(y.ticks(5))
				.bottom(y)
				.strokeStyle(function(d) { return d ? "#eee" : "#000"; })
			  .anchor("left").add(pv.Label)
				.text(y.tickFormat);
			
			this.survbox.add(pv.Label)
				.data([this.y_label])
				.left(-20).bottom(h/2)
				.font("10pt Arial").textAlign("center")
				.textAngle(-Math.PI/2);

			this.survbox.add(pv.Label)
				.data([this.x_label])
				.left(w/2).bottom(-35)
				.font("10pt Arial").textAlign("center");
			
			this.low_censored = this.survbox.add(pv.Dot)
				.shape('cross').strokeStyle('#888').size(8)
				.visible(function(d) { return d.censored; })
				.left(function(d) { return x(d.x); })
				.bottom(function(d) { return y(d.y); });
				
			this.high_censored = this.survbox.add(pv.Dot)
				.shape('cross').strokeStyle('#888').size(8)
				.visible(function(d) { return d.censored; })
				.left(function(d) { return x(d.x); })
				.bottom(function(d) { return y(d.y); });
				
			this.low_line = this.low_censored.add(pv.Line)
				.interpolate("step-after").visible(true)
				.left(function(d) { return x(d.x); })
				.bottom(function(d) { return y(d.y); })
				.strokeStyle(colors[1])
				.lineWidth(2);
				
			this.high_line = this.high_censored.add(pv.Line)
				.interpolate("step-after").visible(true)
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
		
		this.updateHistogram();
		this.updateSurvival();
	}
		
	this.setStatsDiv = function(div) {
		this.stats_div = div;
	}
}

