%rebase base path=path, url=url, scripts='ui_d3,fisheye'

%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title']}}</h2>

<h3>Graph visualization</h3>
<div id="fig"></div>

<div id="sidebar">
<div class="sidebox">
<p><b>Visualization</b><br>
<!--Gene: <input type="text" id="GeneSelect" value="" />-->
</p><br>

<b>Platforms</b><br>
{{platform if 'platform' in locals() else '-'}}
<br><br>

<b>Export</b><br>
<!--<button id="export_table">Table</button>-->
<button id="export_svg">SVG</button>
<br>
</div>
</div>

%if defined('description'):
<h3>Description</h3>
{{description}}
%end

<script type="text/javascript">
var data_file = '{{data}}';

$(document).ready(function() {
	
	var w = 650, h = 500;
	var node_r = 4.5;
	
	$('#export_svg').button().click(function() { export_svg('graph.svg'); });
	
	var fisheye = d3.fisheye().radius(120);

	var force = d3.layout.force()
		.charge(-80).linkDistance(40).size([w, h]);

	var svg = d3.select("#fig").append("svg")
		.attr("width", w).attr("height", h);

	/*svg.append("rect")
		.attr("class", "background")
		.attr("width", width)
		.attr("height", height);*/

	d3.json(data_file, function(data) {
		var n = data.nodes.length;

		force.nodes(data.nodes).links(data.links);

		// Initialize the positions deterministically, for better results.
		data.nodes.forEach(function(d, i) { d.x = d.y = w / n * i; });

		// Run the layout a fixed number of times.
		// The ideal number of times scales with graph complexity.
		force.start();
		for (var i = n*2; i > 0; --i) force.tick();
		force.stop();

		// Center the nodes in the middle.
		var ox = 0, oy = 0;
		data.nodes.forEach(function(d) { ox += d.x, oy += d.y; });
		ox = ox / n - w / 2, oy = oy / n - h / 2;
		data.nodes.forEach(function(d) { d.x -= ox, d.y -= oy; });

		var link = svg.selectAll(".link")
			.data(data.links)
			.enter().append("line")
			.attr("class", "link")
			.attr("x1", function(d) { return d.source.x; })
			.attr("y1", function(d) { return d.source.y; })
			.attr("x2", function(d) { return d.target.x; })
			.attr("y2", function(d) { return d.target.y; })
			.style('stroke', function(d)
				{ return d.type == 1 ? '#999' : '#f88'; })
			.style('stroke-opacity', 0.6)
			.style('stroke-width', function(d) {
				return Math.max(3*d.value, 2); });

		var sel = svg.selectAll(".node")
			.data(data.nodes)
			.enter();

		var node = sel.append("circle")
			.attr("class", "node")
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; })
			.attr("r", node_r)
			.style('stroke', '#fff')
			.style('stroke-width', '1.5px')
			.style("fill", '#0f0');
			//.call(force.drag);
			
		var labels = sel.append('text')
			.attr('text-anchor', 'middle')
			.attr('x', function(d) { return d.x; })
			.attr('y', function(d) { return d.y - 6; })
			.style('font-size', '7pt')
			.text(function(d) { return d.name; });

			
		/*
		var label = svg.append('text')
			.attr('text-anchor', 'middle');
		
		node
			.on('mouseout', function(d) { label.text(''); })
			.on('mouseover', function(d) {
				label.attr('x', d.x).attr('y', d.y - 10).text(d.name);
			});
		*/

			
		svg.on("mousemove", function() {
			fisheye.center(d3.mouse(this));

			node.each(function(d) { d.display = fisheye(d); })
				.attr("cx", function(d) { return d.display.x; })
				.attr("cy", function(d) { return d.display.y; })
				.attr("r", function(d) { return d.display.z * node_r; });
			
			labels.each(function(d) { d.display = fisheye(d); })
				.attr("x", function(d) { return d.display.x; })
				.attr("y", function(d) { return d.display.y - 6; });
				
			link
				.attr("x1", function(d) { return d.source.display.x; })
				.attr("y1", function(d) { return d.source.display.y; })
				.attr("x2", function(d) { return d.target.display.x; })
				.attr("y2", function(d) { return d.target.display.y; });
		});
	});

});
</script>

