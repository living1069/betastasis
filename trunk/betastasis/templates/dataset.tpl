%rebase base path=path, url=url

<h2>{{get('title', path[path.rfind(' > ')+3:])}}</h2>

{{!get('content', '')}}

%subs = [(k, v) for k,v in hierarchy.iteritems() if k.startswith(url + '/') and not '/' in k[len(url)+1:]]

%prev_cat = ''
%for url, val in subs:
	%category = val[2]['category']
	%template = val[1]
	%if prev_cat != category:
		%if prev_cat != '':
</ul>
		%end
		%prev_cat = category
<h3>{{category}}</h3><ul>
	%end
	
	%name = val[0]
	%visname = name[name.rfind(' > ')+3:]
	%if template == 'igv_track':
	<li><a href="http://www.broadinstitute.org/igv/projects/current/igv.php?genome=hg19&sessionURL={{val[2]['data']}}">{{visname}}</a></li>
	%else:
	<li><a href="/{{url}}/">{{visname}}</a></li>
	%end
%end

<div id="sidebar">
<div class="sidebox">
<p><b>Dataset overview</b><br />

%if defined('pubmed_id'):
PubMed ID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/{{pubmed_id}}">{{pubmed_id}}</a><br>
%end

%if defined('geo_series_id'):
GEO series: <a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={{geo_series_id}}">{{geo_series_id}}</a><br>
%end

%if defined('publication_date'):
Published: {{publication_date}}<br>
%end

{{!get('stats', '')}}
</p>
	
%if defined('platform'):
<p><b>Platforms used</b><br>
{{!get('platform', '')}}
</p>
%end

</ul></div></div>

