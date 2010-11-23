// Betastasis math
var bmath = {};

bmath.quartiles = function(p_array)
{
    var array = p_array.slice(0);
    array.sort(function(a,b){return a - b});
    
    var median_idx;
    var median;
    
    var lower = [];
    var lq;

    var upper = [];
    var uq;
    
    if (array.length % 2 > 0)
    {	
	median_idx = Math.floor(array.length / 2);
	median = array[median_idx];

	lower = array.slice(0, median_idx);
	upper = array.slice(median_idx + 1);

	quartile_idx = Math.floor(median_idx / 2);
	lq = (lower[(quartile_idx - 1)] + lower[(quartile_idx)]) / 2;
	uq = (upper[(quartile_idx - 1)] + upper[(quartile_idx)]) / 2;
    }
    else
    {
	median_idx = array.length / 2;
	median = (array[median_idx - 1] + array[median_idx]) / 2;

	lower = array.slice(0, median_idx);
	upper = array.slice(median_idx);
	
	quartile_idx = Math.floor(median_idx / 2);
	lq = (lower[(quartile_idx - 1)] + lower[(quartile_idx)]) / 2;
	uq = (upper[(quartile_idx - 1)] + upper[(quartile_idx)]) / 2;
    }

    var foo = { median: median,
		lq: lq,
		uq: uq };

    return foo;
};

bmath.magn_round = function(d, up)
{	
    var sign = 1;
    if (d < 0)
    {
	d = d * -1;
	sign = -1;
    }
    
    var magn = Math.pow(10, Math.floor(pv.log(d, 10)));
    var foo = d / magn;
    var rounded;
    
    if (up == false)
    {
	rounded = sign * magn*(Math.floor(foo));
    }
    else
    {
	rounded = magn*(Math.ceil(foo));
    }
    
    rounded = sign*rounded;
    return rounded;
};

bmath.cosh = function(x) {
    return 0.5 * (Math.exp(x) + Math.exp(-x));
}
