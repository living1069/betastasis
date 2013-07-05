function Scatterplot(div) {
	var w = 600, h = 600;
	var self = this;
	self.xlabel = 'X';
	self.ylabel = 'Y';
	self.scale = 'natural';
	
	var vis = d3.select(div).append('svg:svg')
		.attr('width', w+50).attr('height', h+50);
	
	var box = vis.append('svg:g')
		.attr('transform', 'translate(10,10)');
	
	this.set_data = function(data) {
		
		if (this.scale != 'natural')
			alert('Only natural scale supported for now.');
		
		self.data = data;
		
		var dx = (d3.max(data.x) - d3.min(data.x)) / 10;
		var dy = (d3.max(data.y) - d3.min(data.y)) / 10;
		
		/*
		var x = d3.fisheye.scale(d3.scale.linear)
			.domain([d3.min(data.x)-dx, d3.max(data.x)+dx]).range([0, w]);
		var y = d3.fisheye.scale(d3.scale.linear)
			.domain([d3.min(data.y)-dy, d3.max(data.y)+dy]).range([h, 0]);
		*/

		var x = d3.scale.linear()
			.domain([d3.min(data.x)-dx, d3.max(data.x)+dx]).range([0, w]);
		var y = d3.scale.linear()
			.domain([d3.min(data.y)-dy, d3.max(data.y)+dy]).range([h, 0]);
		
		box.selectAll('axis').remove();
		var x_axis = box.append('svg:g').classed('axis', true)
			.attr('transform', 'translate(0,' + y(0) + ')');
		x_axis.call(d3.svg.axis().scale(x).orient('bottom'));
	
		var y_axis = box.append('svg:g').classed('axis', true)
			.attr('transform', 'translate(' + x(0) + ',0)');
		y_axis.call(d3.svg.axis().scale(y).orient('left'));
			
		var dots = box.selectAll('.scatter_dots')
			.data(d3.range(data.x.length))
			.enter().append('svg:circle')
				.classed('scatter_dots', true)
				.attr('cx', function (d, i) { return x(data.x[i]); })
				.attr('cy', function (d, i) { return y(data.y[i]); })
				.attr('r', function (d, i) { return data.radius[i]; })
				.style('fill', function (d, i) { return data.color[i]; });
				
		var dot_label = box.append('svg:text')
			.attr('text-anchor', 'middle')
			.style('font-size', '10pt').text('');

		dots
			.on('mouseover', function() {
				var dot = d3.select(this);
				dot_label
					.attr('x', dot.attr('cx'))
					.attr('y', dot.attr('cy') - 10)
					.text(data.label[dot.data()]);
			})
			.on('mouseout', function() {
				dot_label.text('');
			});
		
		/*
		vis.on('mousemove', function() {
			var mouse = d3.mouse(this);
			x.distortion(1).focus(mouse[0]);
			y.distortion(1).focus(mouse[1]);
			
			dots.attr('cx', function (d, i) { return x(data.x[i]); })
				.attr('cy', function (d, i) { return y(data.y[i]); });
			
			x_axis.attr('transform', 'translate(0,' + y(0) + ')')
				.call(d3.svg.axis().scale(x).orient('bottom'));
			y_axis.attr('transform', 'translate(' + x(0) + ',0)')
				.call(d3.svg.axis().scale(y).orient('left'));
		});*/
	}
}

