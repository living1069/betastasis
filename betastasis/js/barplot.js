// Author: Matti Annala <matti.annala@tut.fi>

function Barplot(div) {
    var w = 600;
    var h = 200;
    
    var y_rule_colors = ["#000", "#ccc"];
    var box_color = "lightgray";
    
    var ylabel_width = 50;
	this.diag_labels = true;
	this.scale = 'natural';
	
    this.vis = new pv.Panel()
		.width(w+100).height(h).margin(20).bottom(30)
		.canvas(div);
    
    this.set_data = function(expr, labels) {
		
		// First we discard all samples with missing expression values.
		labels = labels.filter(function(d, i) { return expr[i] != -100; });
		expr = expr.filter(function(d, i) { return expr[i] != -100; });
	
		if (this.scale == 'natural') {
			expr = expr.map(function(d) { return Math.pow(2, d); });
		} else if (this.scale == 'log-10') {
			expr = expr.map(function(d) { return 0.3010 * d; });
		}
		
		// Resolve global max
		var max_expr = -1000;
		for (var s in expr)
			max_expr = (expr[s] > max_expr) ? expr[s] : max_expr;
		
		var plot_min = (this.scale == 'natural') ? 0 : Math.floor(pv.min(expr));
		var plot_max = Math.ceil(1.1 * max_expr);

		// Scale for highlighting.
		// This scale is also used for the boxes with little tweaks.
		var hlx = pv.Scale.ordinal(labels).splitBanded(ylabel_width, w, 1.0);
		var y = pv.Scale.linear(plot_min, plot_max).nice().range(0, h);
		var s = hlx.range().band / 5;
		var bar_width = s * 2;

		// Scale for expression dots in individual bars
		var dot_size = 5;
		var bar_x = pv.Scale.linear(0, 1).range(dot_size, bar_width - dot_size);
		
		/* Add the y-axis rules */
		if (this.vis.y_rules == undefined) {
			this.vis.y_rules = this.vis.add(pv.Rule)
				.left(ylabel_width).width(w-ylabel_width);
			
			this.vis.bottom_rule = this.vis.add(pv.Rule)
				.bottom(0).height(0).left(ylabel_width).width(w-ylabel_width)
				.strokeStyle("#000");
			
			this.vis.y_rules.labels = this.vis.y_rules.anchor("left")
				.add(pv.Label)
		}
		
		this.vis.y_rules.data(y.ticks())
			.bottom(y)
			.strokeStyle("lightgray");

		this.vis.y_rules.labels.text(y.tickFormat);

		if (this.vis.gene_label == undefined) {
			this.vis.gene_label = this.vis.add(pv.Label)
				.bottom(h / 2)
				.font("16pt Arial").textAlign("center").textAngle(-Math.PI/2);
		}
		
		this.vis.gene_label.data([this.gene_name]);
		
		if (this.type_labels == undefined) {
			this.type_labels = this.vis.add(pv.Label);
		}
		
		this.type_labels.data(labels)
			.left(function(d) { return hlx(d) + hlx.range().band / 2 - 7; })
			.textAlign('right');
		
		if (this.diag_labels) {
			if (labels.length > 10) {
				this.type_labels.bottom(-5)
					.textAngle(-Math.PI / 2)
					.left(function(d) {
						return hlx(d) + hlx.range().band / 2 + 7; });
			} else {
				this.type_labels.textAngle(-Math.PI / 4).bottom(-15)
					.left(function(d) {
						return hlx(d) + hlx.range().band / 2 + 10; });
			}
		} else {
			this.type_labels.textAngle(0).bottom(-25);
		}
		
		var xtick_font_size = 12;
		if (labels.length > 40) {
			xtick_font_size = 7;
		} else if (labels.length > 20) {
			xtick_font_size = 10;
		}
		
		this.type_labels.font(xtick_font_size + 'pt Arial');
		


	
		/* Add a panel which will hold the boxes (bars) for each data point */
		if (this.expr_boxes == undefined) {	    
			this.expr_boxes = this.vis.add(pv.Bar)
				.strokeStyle('black').lineWidth(1).antialias(false);
		}
		
		this.type_labels.data(labels);
		
		this.expr_boxes.data(expr)
			.width(bar_width)
			.left(function(d) { return hlx(labels[this.index]) + hlx.range().band / 2 - s; })
			.bottom(y(plot_min)).height(function(d){return y(d) - y(plot_min);})
			.fillStyle(box_color);
		
		// Calculate an estimate for the required X-axis label height.
		var xlabel_lens = labels.map(function(d) { return d.length; });
		var xlabel_width = pv.max(xlabel_lens) * xtick_font_size;
		this.vis.bottom(xlabel_width);
		
		this.vis.render();
    }
}
