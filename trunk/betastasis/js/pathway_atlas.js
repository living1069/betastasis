Array.prototype.abs_max = function() {
	var max = this[0];
	for (var k = 1; k < this.length; k++) {
		if (Math.abs(this[k]) > Math.abs(max))
			max = this[k];
	}
	return max;
}

var annotations_loaded = function(d) {
	
	var annotations = d;
	var svg = document.getElementById('test_svg').getSVGDocument();
	
	var inner_shades = [
		'stop-color:hsl(0,99%,85%); stop-opacity:1',
		'stop-color:hsl(0,66%,85%); stop-opacity:1',
		'stop-color:hsl(0,33%,85%); stop-opacity:1',
		'stop-color:hsl(0,0%,85%); stop-opacity:1',
		'stop-color:hsl(220,33%,85%); stop-opacity:1',
		'stop-color:hsl(220,66%,85%); stop-opacity:1',
		'stop-color:hsl(220,99%,85%); stop-opacity:1'
	];
	
	var outer_shades = [
		'stop-color:hsl(0,99%,70%); stop-opacity:1',
		'stop-color:hsl(0,66%,70%); stop-opacity:1',
		'stop-color:hsl(0,33%,70%); stop-opacity:1',
		'stop-color:hsl(0,0%,70%); stop-opacity:1',
		'stop-color:hsl(220,33%,70%); stop-opacity:1',
		'stop-color:hsl(220,66%,70%); stop-opacity:1',
		'stop-color:hsl(220,99%,70%); stop-opacity:1'
	];
	
	var groups = svg.getElementsByTagName('g');
	for (var k = 0; k < groups.length; k++) {
		var group = groups[k];
		if (group.childNodes.length != 5) continue;
		var text = group.childNodes[3];
		if (text.constructor.toString().indexOf('SVGTextElement') == -1)
			continue;
		
		var symbol = text.textContent;
		var annot = annotations[symbol];
		
		// Color the symbol.
		var level = 0;
		if (annot != undefined && 'levels' in annot) {
			level = annot.levels.abs_max();
		}
		
		var style = group.childNodes[1].style.cssText;
		var match = /(radialGradient[\d-+]+)/.exec(style);
		if (match.length == 2) {
			var rad_grad = svg.getElementById(match[1]);
			var lin_grad = svg.getElementById(rad_grad.href.animVal.slice(1));
			var stop_1 = lin_grad.childNodes[1];
			var stop_2 = lin_grad.childNodes[3];
			var idx = 6 - (level + 3);
			stop_1.setAttributeNS(null, 'style', inner_shades[idx]);
			stop_2.setAttributeNS(null, 'style', outer_shades[idx]);
		}

		if (annot == undefined) continue;

		var desc = '<b>' + annot.name + '</b><br>';
		if ('statements' in annot) {
			desc += '<ul>';
			for (var s = 0; s < annot.statements.length; s++) {
				desc += '<li>' + annot.statements[s] + '</li>';
			}
			desc += '</ul>';
		}
		
		// NOTE: In-loop closure workaround.
		group.onmouseover = (function(val) {
			return function() { tooltip.show(val); }; })(desc);
		group.onmouseout = function() { tooltip.hide(); };
	}
};

var tooltip=function(){
	var id = 'tt';
	var top = 3;
	var left = 3;
	var maxw = 300;
	var speed = 10;
	var timer = 20;
	var endalpha = 95;
	var alpha = 0;
	var tt,t,c,b,h;
	var ie = document.all ? true : false;
	return{
		show:function(v,w){
			if(tt == null){
				tt = document.createElement('div');
				tt.setAttribute('id',id);
				t = document.createElement('div');
				t.setAttribute('id',id + 'top');
				c = document.createElement('div');
				c.setAttribute('id',id + 'cont');
				b = document.createElement('div');
				b.setAttribute('id',id + 'bot');
				tt.appendChild(t);
				tt.appendChild(c);
				tt.appendChild(b);
				document.body.appendChild(tt);
				tt.style.opacity = 0;
				tt.style.filter = 'alpha(opacity=0)';
				var svg = document.getElementById('test_svg').getSVGDocument();
				svg.onmousemove = this.pos;
			}
			tt.style.display = 'block';
			c.innerHTML = v;
			tt.style.width = w ? w + 'px' : 'auto';
			if(!w && ie){
				t.style.display = 'none';
				b.style.display = 'none';
				tt.style.width = tt.offsetWidth;
				t.style.display = 'block';
				b.style.display = 'block';
			}
			if(tt.offsetWidth > maxw){tt.style.width = maxw + 'px'}
			h = parseInt(tt.offsetHeight) + top;
			clearInterval(tt.timer);
			tt.timer = setInterval(function(){tooltip.fade(1)},timer);
		},
		pos:function(e){
			var svg_pos = $('#test_svg').offset();
			var u = svg_pos.top + e.pageY;
			var l = svg_pos.left + e.pageX;
			tt.style.top = (u - h) + 'px';
			tt.style.left = (l + left) + 'px';
		},
		fade:function(d){
			var a = alpha;
			if((a != endalpha && d == 1) || (a != 0 && d == -1)){
				var i = speed;
				if(endalpha - a < speed && d == 1){
					i = endalpha - a;
				}else if(alpha < speed && d == -1){
					i = a;
				}
				alpha = a + (i * d);
				tt.style.opacity = alpha * .01;
				tt.style.filter = 'alpha(opacity=' + alpha + ')';
			}else{
				clearInterval(tt.timer);
				if(d == -1){tt.style.display = 'none'}
			}
		},
		hide:function(){
			clearInterval(tt.timer);
			tt.timer = setInterval(function(){tooltip.fade(-1)},timer);
		}
	};
}();

