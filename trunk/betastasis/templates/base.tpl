<head>
	<title>Project Betastasis</title>
	<meta http-equiv="content-type" content="text/html; charset=UTF-8">
	<script type="text/javascript" src="/js/jquery.js"></script>
	<script type="text/javascript" src="/js/utils.js"></script>
	
%if 'scripts' in locals():
	%for script in scripts.split(','):
		%if script == 'ui':
	<script type="text/javascript" src="/js/protovis.js"></script>
	<script type="text/javascript" src="/js/jquery-ui.js"></script>
	<script type="text/javascript" src="/js/autocomplete.js"></script>
	<script type="text/javascript" src="/js/math_util.js"></script>
	<script type="text/javascript" src="/js/group_selector.js"></script>
		%elif script == 'ui_d3':
	<script type="text/javascript" src="/js/d3.v2.min.js"></script>
	<script type="text/javascript" src="/js/science/science.min.js"></script>
	<script type="text/javascript" src="/js/science/science.stats.min.js"></script>
	<script type="text/javascript" src="/js/jquery-ui.js"></script>
	<script type="text/javascript" src="/js/autocomplete.js"></script>
	<script type="text/javascript" src="/js/math_util.js"></script>
	<script type="text/javascript" src="/js/group_selector.js"></script>

		%elif script == 'ui_circvis':
	<script type="text/javascript" src="/js/autocomplete.js"></script>
	<script type="text/javascript" src="/js/math_util.js"></script>
	<script type="text/javascript" src="/js/group_selector.js"></script>
    <script type="text/javascript" src="/js/visquick-core-1.0/protovis-r3.3.1.js"></script>
    <script type="text/javascript" src="/js/jquery-ui.js"></script>
    <script type="text/javascript" src="/js/visquick-core-1.0/visquick-utils.js"></script>
    <script type="text/javascript" src="/js/circvis/circvis.js"></script>
	<script type="text/javascript" src="/js/cytoband.hg18.json"></script>
	<link rel="stylesheet" type="text/css" href="/style/multiselect.css">

		%elif script == 'chrom_data':
	<script type="text/javascript" src="/js/circvis/data/cytoband.hg19.json"></script>
	<script type="text/javascript" src="/js/circvis/data/chromInfo.json"></script>

		%elif script == 'multiselect':
	<script type="text/javascript" src="/js/jquery.multiselect.min.js"></script>
	<link rel="stylesheet" type="text/css" href="/style/jquery.multiselect.css">

		%else:
	<script type="text/javascript" src="/js/{{script}}.js"></script>
		%end
	%end
%end
	
	<link rel="stylesheet" href="/style/style.css" type="text/css">
	<link rel="stylesheet" type="text/css" href="/style/jquery-ui/jquery-ui.css">
		
</head>

<!-- Test for browser SVG support -->
<script type="text/javascript">
if (!document.implementation.hasFeature("http://www.w3.org/TR/SVG11/feature#BasicStructure", "1.1")) {
	alert('Your browser lacks support for the visualizations on this website. Supported browsers include Google Chrome, Mozilla Firefox and Internet Explorer 9 or later.');
}
</script>

<body>
<div id="header">
	<h1 id="logo"><a href="/">Project Betastasis</a></h1>
</div>
	
<ul class="xbreadcrumbs" id="breadcrumbs">
	<li><a href="/">Home</a></li>
%if url:
	%levels = zip(path.split(' > '), url.split('/'))
	%curr_level = '/'
	%for name, level in levels:
		<li><a href="{{curr_level+level}}/">{{name}}</a></li>
		%curr_level += level + '/'
	%end
%end
</ul>

<div id="main">
%include
</div>

</body></html>


<!-- GOOGLE ANALYTICS -->
<script type="text/javascript">
  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-29086022-1']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();
</script>
