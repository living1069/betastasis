<head>
	<meta http-equiv="content-type" content="text/html; charset=UTF-8">
	<script type="text/javascript" src="/js/jquery.js"></script>
	<script type="text/javascript" src="/js/xbreadcrumbs.js"></script>
	<script type="text/javascript" src="/js/utils.js"></script>
	
%if 'scripts' in locals():
	%for script in scripts:
		%if script == 'ui':
	<script type="text/javascript" src="/js/protovis.js"></script>
	<script type="text/javascript" src="/js/jquery-ui.js"></script>
	<script type="text/javascript" src="/js/autocomplete.js"></script>
	<script type="text/javascript" src="/js/math_util.js"></script>
	<script type="text/javascript" src="/js/sprintf.js"></script>
	<script type="text/javascript" src="/js/group_selector.js"></script>
		%else:
	<script type="text/javascript" src="/js/{{script}}.js"></script>
		%end
	%end
%end
	
	<link rel="stylesheet" href="/style/style.css" type="text/css">
	<link rel="stylesheet" href="/style/xbreadcrumbs.css" />
	<link rel="stylesheet" type="text/css" href="/style/jquery-ui/jquery-ui.css">
		
	
	<title>Project Betastasis</title>
	
	<script type="text/javascript">
		$(document).ready(function() { $('#breadcrumbs').xBreadcrumbs(); });
	</script>
</head>
<body>

<div id="header">
	<h1 id="logo"><a href="/">Project Betastasis</a></h1>
</div>
	
<ul class="xbreadcrumbs" id="breadcrumbs">
	<li><a href="/">Home</a></li>
%levels = nav_tree.pop('levels')
%node = nav_tree
%curr_level = '/'
%for level in levels:
	%v = node[level]
	<li><a href="{{curr_level+level}}/">{{v['index'] if type(v) is dict else v}}</a>\\
	%if len(node) <= 2:
		</li>
	%else:
<ul>
		%keys = node.keys()
		%keys.sort()
		%for k in keys:
			%v = node[k]
			%if k == 'index': continue
			%name = v['index'] if type(v) is dict else v
		<li><a href="{{curr_level+k}}/">{{name}}</a></li>
		%end
	</ul></li>
	%end
	
	%curr_level += level + '/'
	%node = node[level]
%end
</ul>

%include

</body></html>
