function ExprBoxplot(div, opts) {
    // Plot variables
    var w = 550;
    var h = 200;
    
    // Colors
    var y_rule_colors = ["#000", "#ccc"];
    var box_color = "lightgray";
    var highlight_color = pv.color("black").alpha(0.1);
    
    // Gene name label width
    var label_width = 40;

    // Precision of values shown in tooltips
    var tooltip_prec = 3;
    
    // Root panel
    var vis = new pv.Panel()
	.width(w)
	.height(h)
	.margin(20)
	.canvas(div);
 
    // Parse options
    this.show_expr_dots = false;
    if (opts != undefined) {
	if (opts.show_dots != undefined) {
	    this.show_expr_dots = opts.show_dots;
	}
    }
    
    this.vis = vis;
    this.visible_panel = new pv.Panel()
	.visible(true);

    // Gene name label on the y-axis
    this.vis.gene_label = undefined;
    
    // Y axis rules
    this.vis.y_rules = undefined;
    
    // Labels
    this.type_labels = undefined;
    
    // Boxes
    this.expr_boxes = undefined;
    this.expr_box_visible = true;
    
    // Expression data dots
    this.expr_dot_panel = undefined;
    this.expr_dot_dots = undefined;
    this.expr_dot_xcoord = undefined;
    this.expr_dot_visible = false;
    
    // Mouseover
    this.expr_highlights = undefined;
    
    this.set_data = function(param_data) {
	var subtype_expr = param_data.contents.expr_data;

	// Calculate random coordinates for expression dots
	for (var i = 0; i < subtype_expr.length; i = i + 1) {
	    var expr_dots = [];
	    var expr = subtype_expr[i].expr_values;
	    for (var j = 0; j < expr.length; j = j + 1) {		
		var relative_pos = (expr[j] - subtype_expr[i].median) / (subtype_expr[i].max - subtype_expr[i].min) * 1.5;		
		expr_dots[j] = { value: expr[j],
				 x: ((Math.random()-0.5) * Math.abs(2-bmath.cosh(relative_pos*Math.PI/2))) + 0.5
			       };
	    }
	    subtype_expr[i].expr_dots = expr_dots;
	}
	
	// Resolve global min and max
	var min_expr = pv.min(subtype_expr.map(function (d) {
	    return d.min;
	}));

	var max_expr = pv.max(subtype_expr.map(function (d) {
	    return d.max;
	}));

	var expr_range = max_expr - min_expr;
	var plot_min = Math.floor(min_expr - 0.1 * expr_range);
	var plot_max = Math.ceil(max_expr + 0.1 * expr_range);

	// Scale for highlighting.
	// This scale is also used for the boxes with little tweaks.
	var hlx = pv.Scale.ordinal(subtype_expr, function (g) {return g.subtype}).splitBanded(label_width, w, 1.0);
	var y = pv.Scale.linear(plot_min, plot_max).nice().range(0, h);
	var s = hlx.range().band / 5;
	var bar_width = s * 2;

	// Scale for expression dots in individual bars
	var dot_size = 5;
	var bar_x = pv.Scale.linear(0, 1).range(dot_size, bar_width - dot_size);
	
	/* Add the y-axis rules */
	if (this.vis.y_rules === undefined) {
	    this.vis.y_rules = this.vis.add(pv.Rule)
		.left(label_width)
		.width(w-label_width)
		
	    this.vis.bottom_rule = this.vis.add(pv.Rule)
		.bottom(0)
		.height(0)
		.left(label_width)
		.width(w-label_width)
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
		.font("16pt Arial")
		.textAlign("center")
		.textAngle(-Math.PI / 2);
	}
	
	this.vis.gene_label.data([param_data.contents.gene_name]);
	
	if (this.type_labels === undefined) {
	    this.type_labels = this.vis.add(pv.Label)
		.left(function (d){return hlx(d.subtype) + hlx.range().band / 2})
		.textAlign("center")
		.bottom(-20)
		.font("12pt Arial");
	}
	
	/* Add a panel which will hold the boxes (bars) for each data point */
	if (this.expr_boxes === undefined) {	    
	    this.expr_boxes = this.vis.add(pv.Panel)
		.left(function (d){return hlx(d.subtype) + hlx.range().band / 2 - s})
		.visible(this.expr_box_visible)
		.width(bar_width);
	    
	    // Range line of the box
	    this.expr_boxes.range_line = this.expr_boxes.add(pv.Rule)
		.left(s);
	    
	    // Min and max indicators
	    this.expr_boxes.minmax = this.expr_boxes.add(pv.Rule)
		.left(s / 2)
		.width(s);
	    
	    // Upper/lower quartile ranges
	    this.expr_boxes.quartiles = this.expr_boxes.add(pv.Bar)
		.strokeStyle("black")
		.lineWidth(1)
		.antialias(false);
	    
	    // Median line
	    this.expr_boxes.median = this.expr_boxes.add(pv.Rule);
	}
	
	this.type_labels.data(subtype_expr)
	    .text(function (d){return d.subtype});
	this.expr_boxes.data(subtype_expr);	
	this.expr_boxes.range_line
	    .bottom(function (d){return y(d.min)})
	    .height(function (d){return y(d.max) - y(d.min)});
	this.expr_boxes.minmax
	    .data(function (d){return [d.min, d.max]})
	    .bottom(y);
	this.expr_boxes.quartiles
	    .bottom(function (d){return y(d.lq)})
	    .height(function (d){return y(d.uq) - y(d.lq)})
	    .fillStyle(function (d){return box_color});
	this.expr_boxes.median.bottom(function (d){return y(d.median)});
	
	if (this.expr_highlights === undefined)
	{
	    this.expr_highlights = this.vis.add(pv.Panel)
		.left(function (d){return hlx(d.subtype)})
		.def("i", -1)
		.fillStyle(function() {return this.i() == this.index ? highlight_color : pv.color("white").alpha(0.01)})
		.event("mouseover", pv.Behavior.tipsy({gravity: "n", fade: true, html: true}))
		.event("mouseout", function() {return this.i(-1)})
		.text(this.tooltip_text)
		.text(function(d) {
		    return "<table><tr><td>" + 
			"<label class=\"tooltip_label\">max:</label></td><td><label class=\"tooltip_field\">" + d.max.toPrecision(tooltip_prec) + "</label></td></tr><tr><td>" +
			"<label class=\"tooltip_label\">uq:</label></td><td><label class=\"tooltip_field\">" + d.uq.toPrecision(tooltip_prec) + "</label></tr><tr><td>" +
			"<label class=\"tooltip_label\">med:</label></td><td><label class=\"tooltip_field\">" + d.median.toPrecision(tooltip_prec) + "</label></td></tr><tr><td>" +
			"<label class=\"tooltip_label\">lq:</label></td><td><label class=\"tooltip_field\">" + d.lq.toPrecision(tooltip_prec) + "</label></td></tr><tr><td>" +
			"<label class=\"tooltip_label\">min:</label></td><td><label class=\"tooltip_field\">" + d.min.toPrecision(tooltip_prec) +
			"</td></tr></table>"})
		.width(hlx.range().band);	    
	}

	this.expr_highlights.data(subtype_expr);
	
	if (this.expr_dot_panel === undefined) {
	    this.expr_dot_panel = this.vis.add(pv.Panel)
		.left(function (d){return hlx(d.subtype) + hlx.range().band / 2 - s})
		.visible(this.expr_dot_visible)
		.width(bar_width);
	    
	    this.expr_dot_dots = this.expr_dot_panel.add(pv.Dot)
		.shape('circle')
		.size(2)
		.lineWidth(1)
		.strokeStyle( pv.color('black').alpha(.5) )
		.fillStyle(function() {return this.strokeStyle().alpha(1)});
	}
	
	this.expr_dot_panel
	    .data(subtype_expr);
	
	this.expr_dot_dots
	    .data(function(d) { return d.expr_dots })
	    .bottom(function(d) { return y(d.value) })
	    .left(function(d) { return bar_x(d.x) });

	this.vis.render();
    }

    this.set_box_visibility = function(visible_bit) {
	if (this.expr_box_visible !== visible_bit) {
	    this.expr_box_visible = visible_bit;
	    this.expr_boxes.visible(visible_bit);
	    
	    this.vis.render();
	}
    }
    
    this.set_dots_visibility = function(visible_bit) {
	if (this.expr_dot_visible !== visible_bit) {
	    this.expr_dot_visible = visible_bit;
	    this.expr_dot_panel.visible(visible_bit);
	    
	    this.vis.render();
	}
    }
}
