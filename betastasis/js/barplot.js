// Author: Matti Annala <matti.annala@tut.fi>
// Based on gene expression boxplot visualization by Kalle Leinonen.

function Barplot(div) {
    var w = 550;
    var h = 200;
    
    var y_rule_colors = ["#000", "#ccc"];
    var box_color = "lightgray";
    
    var label_width = 60;
	this.diag_labels = true;
	this.scale = 'natural';

    this.vis = new pv.Panel()
		.width(w+100).height(h).margin(20).bottom(30)
		.canvas(div);
    
    this.set_data = function(expr, samples) {
	
		if (this.scale == 'natural') {
			expr = expr.map(function(d) { return Math.pow(2, d); });
		}
		
		// Resolve global max
		var max_expr = -1000;
		for (var s in expr)
			max_expr = (expr[s] > max_expr) ? expr[s] : max_expr;
		
		var plot_max = Math.ceil(1.1 * max_expr);

		// Scale for highlighting.
		// This scale is also used for the boxes with little tweaks.
		var hlx = pv.Scale.ordinal(samples).splitBanded(label_width, w, 1.0);
		var y = pv.Scale.linear(0, plot_max).nice().range(0, h);
		var s = hlx.range().band / 5;
		var bar_width = s * 2;

		// Scale for expression dots in individual bars
		var dot_size = 5;
		var bar_x = pv.Scale.linear(0, 1).range(dot_size, bar_width - dot_size);
		
		/* Add the y-axis rules */
		if (this.vis.y_rules === undefined) {
			this.vis.y_rules = this.vis.add(pv.Rule)
				.left(label_width).width(w-label_width);
			
			this.vis.bottom_rule = this.vis.add(pv.Rule)
				.bottom(0).height(0).left(label_width).width(w-label_width)
				.strokeStyle("#000");
			
			this.vis.y_rules.labels = this.vis.y_rules.anchor("left")
				.add(pv.Label)
		}
		
		this.vis.y_rules.data(y.ticks())
			.bottom(y)
			.strokeStyle("lightgray");

		this.vis.y_rules.labels.text(y.tickFormat);

		if (this.vis.gene_label === undefined) {
			this.vis.gene_label = this.vis.add(pv.Label)
				.bottom(h / 2)
				.font("16pt Arial").textAlign("center").textAngle(-Math.PI/2);
		}
		
		this.vis.gene_label.data([this.gene_name]);
		
		if (this.type_labels === undefined) {
			this.type_labels = this.vis.add(pv.Label)
				.bottom(-20)
				.font("12pt Arial");
		}
		
		if (this.diag_labels) {
			this.type_labels.textAngle(Math.PI / 4);
		} else {
			this.type_labels.textAngle(0).bottom(-25);
		}
	
		/* Add a panel which will hold the boxes (bars) for each data point */
		if (this.expr_boxes === undefined) {	    
			this.expr_boxes = this.vis.add(pv.Bar)
				.strokeStyle("black").lineWidth(1).antialias(false);
		}
		
		this.type_labels.data(samples)
			.left(function(d) { return hlx(d) + hlx.range().band / 2 - 7; });
		this.expr_boxes.data(expr)
			.width(bar_width)
			.left(function(d) { return hlx(samples[this.index]) + hlx.range().band / 2 - s; })
			.bottom(y(0)).height(function(d){return y(d) - y(0);})
			.fillStyle(box_color);
		
		if (this.diag_labels)
			this.vis.bottom(250);
		
		this.vis.render();
    }
}
