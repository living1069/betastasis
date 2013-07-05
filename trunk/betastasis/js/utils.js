// http://www.netlobo.com/url_query_string_javascript.html
function gup(name) {
  name = name.replace(/[\[]/,"\\\[").replace(/[\]]/,"\\\]");
  var regexS = "[\\?&]"+name+"=([^&#]*)";
  var regex = new RegExp( regexS );
  var results = regex.exec( window.location.href );
  if( results == null )
    return "";
  else
    return results[1];
}

var download_data_uri = function(data, filename, type) {
	var data_url = 'data:' + type + ';base64,' +
		btoa(unescape(encodeURIComponent(data)));
		
	var a = document.createElement('a');
	a.href = data_url;
	a.download = filename;
	a.type = type;
	a.target = '_new';
	a.click();
}

var export_svg = function(filename) {
	$('svg').attr({version: '1.1', xmlns:"http://www.w3.org/2000/svg"});
	download_data_uri($('svg').parent().html(), filename, 'image/svg+xml');
}

