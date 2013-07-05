%rebase base path=path, url=url

<h2>{{get('title', path)}}</h2>

{{get('content', '')}}

<h3>Available datasets</h3><ul>
%subs = [(k,v[0]) for k,v in hierarchy.iteritems() if k.startswith(url + '/') and not '/' in k[len(url)+1:]]
%for url, name in subs:
	%if not 'hidden' in hierarchy[url][2]:
	<li><a href="/{{url}}/">{{name[name.rfind(' > ')+3:]}}</a></li>
	%end
%end

%include news

