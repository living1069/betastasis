function Scatterplot(div) {
	var w = 400, h = 400;
	var self = this;
	self.xlabel = 'X';
	self.ylabel = 'Y';
	self.pearson = 0;
	self.scale = 'log2';
	
	this.vis = new pv.Panel()
		.width(w+100).height(h+50)
		.canvas(div);

	this.scatterbox = this.vis.add(pv.Panel)
		.width(w).height(h)
		.left(50).top(5);
	
	self.vis.render();
	
	this.set_data = function(data) {
		var xval, yval;
		
		if (this.scale == 'natural') {
			xval = data.xval.map(function(d) { return Math.pow(2, d); });
			yval = data.yval.map(function(d) { return Math.pow(2, d); });
		} else if (this.scale == 'log2') {
			xval = data.xval;
			yval = data.yval;
		} else if (this.scale == 'log10') {
			xval = data.xval.map(function(d) { return 0.3010 * d; });
			yval = data.yval.map(function(d) { return 0.3010 * d; });
		} else {
			alert('Invalid scale.');
		}
		
		self.xval = xval;
		self.yval = yval;
		
		var dx = (pv.max(xval) - pv.min(xval)) / 10;
		var dy = (pv.max(yval) - pv.min(yval)) / 10;
		var x = pv.Scale.linear(pv.min(xval)-dx, pv.max(xval)+dx).range(0,w);
		var y = pv.Scale.linear(pv.min(yval)-dy, pv.max(yval)+dy).range(0,h);
		
		self.x = x;
		self.y = y;
		
		var xmean = pv.mean(xval);
		var ymean = pv.mean(yval);
		var x_stdev = pv.deviation(xval);
		var y_stdev = pv.deviation(yval);
		
		var pearson = 0;
		for (var k = 0; k < xval.length; k++) {
			pearson += (xval[k] - xmean) * (yval[k] - ymean);
		}
		pearson /= x_stdev * y_stdev * (xval.length - 1);
		self.pearson = Math.round(pearson * 1000) / 1000;
		
		if (this.scatterbox_xticks == undefined) {
			
			this.scatterbox_xticks = this.scatterbox.add(pv.Rule);
			this.scatterbox_xticklabels = this.scatterbox_xticks
			  .anchor("bottom").add(pv.Label);
			this.scatterbox_yticks = this.scatterbox.add(pv.Rule);
			this.scatterbox_yticklabels = this.scatterbox_yticks
			  .anchor("left").add(pv.Label);
			
			this.scatterbox.add(pv.Rule).data([0])
				.strokeStyle('#000').left(0);
			this.scatterbox.add(pv.Rule).data([0])
				.strokeStyle('#000').bottom(0);

			this.scatterbox_ylabel = this.scatterbox.add(pv.Label)
				.left(-35).bottom(h/2)
				.font("10pt Arial").textAlign("center")
				.textAngle(-Math.PI/2);

			this.scatterbox_xlabel = this.scatterbox.add(pv.Label)
				.left(w/2).bottom(-35)
				.font("10pt Arial").textAlign("center");
			
			this.scatterbox_dots = this.scatterbox.add(pv.Dot);
			
			this.scatterbox_fit = this.scatterbox.add(pv.Line);
		}
		
		this.scatterbox_fit.visible(false);
		
		this.scatterbox_dots.data(xval)
			.left(function(d) { return x(xval[this.index]); })
			.bottom(function(d) { return y(yval[this.index]); })
			.size(2).fillStyle('#000')
			.title(function(d) { return data.sample_ids[this.index]; });

		this.scatterbox_xticks.data(x.ticks())
			.left(x).strokeStyle("#eee");
		this.scatterbox_xticklabels.text(x.tickFormat);

		this.scatterbox_yticks.data(y.ticks())
			.bottom(y).strokeStyle("#eee");
		this.scatterbox_yticklabels.text(y.tickFormat);
		
		this.scatterbox_xlabel.data([self.xlabel]);
		this.scatterbox_ylabel.data([self.ylabel]);
		
		self.vis.render();
	}
	
	this.fit_least_squares = function() {
		var xval = self.xval, yval = self.yval;
		
		var xmin = pv.min(xval);
		var ymin = pv.min(yval);
		
		var xmean = pv.mean(xval);
		var ymean = pv.mean(yval);
		
		var sum_x_sq = xval.reduce(
			function(a, b) { return a + b*b; }, 0);
		var sum_xy = xval.reduce(function(a, b, idx) {
			return a + xval[idx] * yval[idx]; }, 0);
		
		var a = (ymean * sum_x_sq - xmean * sum_xy) /
			(sum_x_sq - xval.length * xmean * xmean);
		var b = (sum_xy - xval.length * xmean * ymean) /
			(sum_x_sq - xval.length * xmean * xmean);
		
		var xpos = [pv.min(xval), pv.max(xval)];
		var ypos = [a+b*xpos[0], a+b*xpos[1]];
		
		this.scatterbox_fit.data([[xpos[0], ypos[0]], [xpos[1], ypos[1]]])
			.left(function(d) { return self.x(d[0]); })
			.bottom(function(d) { return self.y(d[1]); })
			.lineWidth(3).strokeStyle('f00').visible(true);
		this.scatterbox_fit.render();
		
		return [a, b];
	}
}

