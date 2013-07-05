%rebase base path=path, url=url

<h2>Select a dataset</h2>

%diseases = [(k,v[0]) for k,v in hierarchy.iteritems() if not '/' in k]
%for disease_url, disease_name in diseases:
<h3><a href="/{{disease_url}}/">{{disease_name}}</a></h3><ul>
	%datasets = [(k,v[0]) for k,v in hierarchy.iteritems() if k.startswith(disease_url + '/') and not '/' in k[len(disease_url)+1:]]
	%for url, name in datasets:
		%if not 'hidden' in hierarchy[url][2]:
	<li><a href="/{{url}}/">{{name[name.rfind(' > ')+3:]}}</a></li>
		%end
	%end
</ul>
%end

<!--

<h3><a href="/tissues/">Tissues</a></h3><ul>
	<li><a href="/tissues/affymetrix/">Affymetrix tissue dataset</a></li>
</ul>

<h3>Platforms</h3><ul>
	<li><a href="/microarrays">Microarrays</a></li>
</ul>

<h3>Tools</h3><ul>
	<li><a href="/tools/seqviewer/">DNA sequence viewer</a></li>
</ul>
-->
		
%include news

<div class="sidebox">
	<p><b>Contact us</b><br>
	<a id="matti_email" href="mailto:first.last@tut.fi">Matti Annala</a><br>
	</p>
</div>

<!-- OBFUSCATE EMAIL -->
<script type="text/javascript">
document.getElementById('matti_email').href = 'mailto:matti.annala' + '@tut.fi';
</script>

