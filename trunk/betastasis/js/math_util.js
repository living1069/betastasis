
var bmath = {};

bmath.quartiles = function(p_array) {
    var array = p_array.slice(0);
    array.sort(function(a,b){return a - b});
    
    var median_idx;
    var median;
    
    var lower = [];
    var lq;

    var upper = [];
    var uq;
    
    if (array.length % 2 > 0) {	
		median_idx = Math.floor(array.length / 2);
		median = array[median_idx];

		lower = array.slice(0, median_idx);
		upper = array.slice(median_idx + 1);

		quartile_idx = Math.floor(median_idx / 2);
		lq = (lower[(quartile_idx - 1)] + lower[(quartile_idx)]) / 2;
		uq = (upper[(quartile_idx - 1)] + upper[(quartile_idx)]) / 2;
    } else {
		median_idx = array.length / 2;
		median = (array[median_idx - 1] + array[median_idx]) / 2;

		lower = array.slice(0, median_idx);
		upper = array.slice(median_idx);
		
		quartile_idx = Math.floor(median_idx / 2);
		lq = (lower[(quartile_idx - 1)] + lower[(quartile_idx)]) / 2;
		uq = (upper[(quartile_idx - 1)] + upper[(quartile_idx)]) / 2;
    }

    return { median: median, lq: lq, uq: uq };
};

/*
Math.magn_round = function(d, up) {
    var sign = 1;
    if (d < 0) {
	d = d * -1;
	sign = -1;
    }
    
    var magn = Math.pow(10, Math.floor(pv.log(d, 10)));
    var foo = d / magn;
    var rounded;
    
    if (up == false) {
	rounded = sign * magn*(Math.floor(foo));
    } else {
	rounded = magn*(Math.ceil(foo));
    }
    
    rounded = sign*rounded;
    return rounded;
};
*/

bmath.cosh = function(x) {
    return 0.5 * (Math.exp(x) + Math.exp(-x));
}





/*
Array.prototype.iterateSorted = function(obj, sort, worker) {
	var keys = [];
    for (var key in obj) {
        if (this.hasOwnProperty[key])
            keys.push(key);
    }
    keys.sort(sort);

    for (var i = 0; i < keys.length; i++) {
        worker(obj[keys[i]], keys[i]);
    }
}
*/

Array.prototype.create = function(type, data) {
	return this.selectAll('.foobar').data(data).enter().append(type);
}

/*
d3.selection.prototype.align_center = function() {
	return this.attr('text-anchor', 'middle')
		.attr('alignment-baseline', 'central');
}
*/

Array.prototype.natSort = function() {
  return this.sort(function(a, b) {
    aa = a.split(/(\d+)/);
    bb = b.split(/(\d+)/);

    for(var x = 0; x < Math.max(aa.length, bb.length); x++) {
      if(aa[x] != bb[x]) {
        var cmp1 = (isNaN(parseInt(aa[x],10)))? aa[x] : parseInt(aa[x],10);
        var cmp2 = (isNaN(parseInt(bb[x],10)))? bb[x] : parseInt(bb[x],10);
        if(cmp1 == undefined || cmp2 == undefined)
          return aa.length - bb.length;
        else
          return (cmp1 < cmp2) ? -1 : 1;
      }
    }
    return 0;
  });
}

// http://stackoverflow.com/questions/3730510/javascript-sort-array-and-return-an-array-of-indicies-that-indicates-the-positi
Array.prototype.sortWithIndices = function() {
  for (var i = 0; i < this.length; i++) {
    this[i] = [this[i], i];
  }
  this.sort(function(left, right) {
    return left[0] < right[0] ? -1 : 1;
  });
  var indices = [];
  for (var j = 0; j < this.length; j++) {
    indices.push(this[j][1]);
    this[j] = this[j][0];
  }
  return indices;
}


