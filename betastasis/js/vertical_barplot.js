// Author: Matti Annala <matti.annala@tut.fi>

function VerticalBarplot(div) {
    var w = 500, barh = 10;
    
    var y_rule_colors = ["#000", "#ccc"];
    var box_color = "lightgray";
    
    var label_width = 20;
	
	page = self;
	self = this;
	self.limits = [];
	self.xlabel = '';

    this.vis = new pv.Panel()
		.top(20)
		.width(w+label_width+50)
		.canvas(div);
	
	this.bar_panel = this.vis.add(pv.Panel)
		.width(w).left(label_width);
	
	this.x_rules = this.bar_panel.add(pv.Rule)
		.strokeStyle('lightgray');
	this.x_rule_labels = this.x_rules.anchor('bottom').add(pv.Label);
	this.x_zero_rule = this.bar_panel.add(pv.Rule).data([0]);
	
	this.x_label = this.bar_panel.add(pv.Label)
		.bottom(-30).textAlign('center')
		.left(w/2);
	
	this.labels = this.bar_panel.add(pv.Label)
		.text(function(d) { return d; })
		.font('10px sans-serif')
		.left(0).top(function(d) { return this.index * (barh + 5) + 16; });
	
	this.bars = this.bar_panel.add(pv.Bar)
		.strokeStyle("black").lineWidth(1).antialias(false);
    
    this.set_data = function(data) {
		var values = data.values;
		var labels = data.labels;
		
		if (self.limits.length == 0) {
			self.limits = [pv.min(values), pv.max(values)];
			self.limits = [
				self.limits[0] - 0.1 * (self.limits[1] - self.limits[0]),
				self.limits[1] + 0.1 * (self.limits[1] - self.limits[0])
			];
		}
		
		var x = pv.Scale.linear(self.limits[0], self.limits[1]).range(0, w);
		var h = (barh + 5) * values.length + 5;
		
		this.vis.height(h + 50);
		this.bar_panel.height(h);
		
		this.x_zero_rule.left(x(0));
		this.x_rules.data(x.ticks())
			.left(function(d) { return x(d); });
		this.x_rule_labels.text(x.tickFormat);
		
		this.labels.data(labels)
			.left(function(d) { 
				return x(values[this.index]) +
					(values[this.index] < 0 ? -3 : 3);
			})
			.textAlign(function(d) {
				return values[this.index] < 0 ? 'right' : 'left'; });
		
		this.bars.data(values)
			.left(function(d) { return d < 0 ? x(d) : x(0); })
			.width(function(d) { return Math.abs(x(0) - x(d)); })
			.top(function(d) { return this.index * (barh + 5) + 5; })
			.height(barh);
		
		if ('urls' in data) {
			this.bars.event('click', function(d) {
				window.open(data.urls[this.index], '_blank');
			});
			this.bars.cursor('pointer');
		}
		
		this.x_label.data([self.xlabel]);
		
		this.vis.render();
    }
}
