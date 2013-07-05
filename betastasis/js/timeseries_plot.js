function TimeseriesPlot(div) {
	var w = 550, h = 300;
	var self = this;
	self.xlabel = 'Time (h)';
	self.ylabel = 'Log-2 ratio vs untreated LNCaP';
	
	var vis = d3.select(div).append('svg:svg')
		.attr('width', w+50).attr('height', h+60);
	
	var box = vis.append('svg:g')
		.attr('transform', 'translate(50,10)');
	
	box.append('svg:line').attr('y2', h).style('stroke', 'black');
	box.append('svg:text')
		.attr('x', -40).attr('y', h/2).text(self.ylabel)
		.attr('text-anchor', 'middle')
		.attr('transform', 'rotate(270,-40,'+(h/2)+')');

	box.append('svg:line').attr('x2', w).attr('y1', h).attr('y2', h)
		.style('stroke', 'black');
	box.append('svg:text')
		.attr('x', w/2).attr('y', h+35).text(self.xlabel)
		.attr('text-anchor', 'middle');
	
	var colors = d3.scale.category10();
	
	this.set_data = function(data) {
		
		var x = d3.scale.linear().range([0,w]);
		var y = d3.scale.linear().range([h,0]);

		x.domain([0, self.max_time]);
		y.domain([Math.min(self.min_expr, -0.2), Math.max(self.max_expr, 0.2)]);
				
		box.selectAll('.axismark').remove();
		
		box.create('svg:text', x.ticks(5)).attr('class', 'axismark')
			.attr('x', x).attr('y', h+18)
			.text(function(d) { return d.toFixed(0); })
			.attr('text-anchor', 'middle');
		
		box.create('svg:line', x.ticks(5)).attr('class', 'axismark')
			.attr('x1', x).attr('x2', x)
			.attr('y1', h).attr('y2', h+5).style('stroke', 'black');
		
		box.create('svg:text', y.ticks(5)).attr('class', 'axismark')
			.attr('x', -18).attr('y', y).attr('dy', 4)
			.text(function(d) { return d.toFixed(1); })
			.attr('text-anchor', 'middle');
			  
		box.create('svg:line', y.ticks(5)).attr('class', 'axismark')
			.attr('x1', -5).attr('x2', 0).attr('y1', y).attr('y2', y)
			.style('stroke', 'black');

		
		// Render a legend.
		box.create('svg:circle', d3.keys(data)).attr('class', 'axismark')
			.attr('r', 4).attr('cx', w - 110)
			.attr('cy', function(d, i) { return 10 + i*20; })
			.style('stroke', 'none')
			.style('fill', function(d, i) { return colors(i); });
		
		box.create('svg:text', d3.keys(data)).attr('class', 'axismark')
			.attr('x', w - 100).attr('text-anchor', 'start')
			.attr('alignment-baseline', 'central')
			.attr('y', function(d, i) { return 10 + i*20; }).text(String);
		
		
		
		
		box.selectAll('.plot').remove();
		
		
		
		var cond_num = 0;
		for (var cond_name in data) {
			var cond = data[cond_name];
			
			var times = d3.keys(cond);
			times = times.map(parseFloat);
			times.sort(d3.ascending);
			
			var means = [];
			var std_of_means = [];
			
			for (var i = 0; i < times.length; i++) {
				var expr = cond[times[i]];
				means.push(d3.mean(expr));
				std_of_means.push(Math.sqrt(
					science.stats.variance(expr) / expr.length +
					self.control_variance_per_N));
			}
			
			var area = d3.svg.area()
				.x(function(d, i) { return x(times[i]); })
				.y0(function(d, i) { return y(d - std_of_means[i]); })
				.y1(function(d, i) { return y(d + std_of_means[i]); });
			
			box.append('svg:path').attr('class', 'plot')
				.attr('d', area(means))
				.style('fill', colors(cond_num))
				.style('fill-opacity', 0.2);
			
			var line = d3.svg.line()
				.x(function(d, i) { return x(times[i]); })
				.y(function(d) { return y(d); });
			
			box.append('svg:path').attr('class', 'plot').attr('d', line(means))
				.style('stroke', colors(cond_num))
				.style('stroke-width', 2).style('fill', 'none');
			
		
			cond_num += 1;
		}
	}
}

