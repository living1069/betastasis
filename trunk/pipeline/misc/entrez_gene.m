function [] = entrez_gene(gene_names, organism)

organism_name = 'Homo sapiens';
if nargin == 2
	organism_name = organism;
end

if iscellstr(gene_names)
	for k = 1:length(gene_names)
		summarize_gene(gene_names{k}, organism_name);
		if k ~= length(gene_names)
			fprintf(1, ['==============================================' ...	
			            '===============================\n']);
			pause(3);
		end
	end
elseif ischar(gene_names)
	summarize_gene(gene_names, organism_name);
else
	error 'Gene names parameter has an invalid type.';
end

return;




function [] = summarize_gene(gene_name, organism_name)

entrez_base_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils';
database = 'gene';

search_results = urlread(sprintf('%s/esearch.fcgi?db=%s&term=%s+AND+(%s[Organism])', ...
	entrez_base_url, database, gene_name, escape_string(organism_name)));

ids = regexp(search_results, '<Id>(\w+)</Id>', 'tokens');
if length(ids) == 0
	fprintf(1, 'No matches found.\n');
	return;
end

best_match = ids{1}; best_match = best_match{1};

%info = urlread(sprintf('%s/efetch.fcgi?db=%s&retmode=xml&id=%s', ...
%	entrez_base_url, database, best_match));

summary_query = sprintf('%s/esummary.fcgi?db=%s&id=%s', entrez_base_url, ...
	database, best_match);

summary = urlread(summary_query);

tokens = regexp(summary, '<Item Name="Name".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, 'GENE NAME:\t%s', tokens{1});

tokens = regexp(summary, '<Item Name="Description".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, ' - %s\n', tokens{1});

tokens = regexp(summary, '<Item Name="OtherAliases".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, 'ALIASES:\t%s\n', tokens{1});

tokens = regexp(summary, ...
	'<Item Name="OtherDesignations".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
descs = strread(tokens{1}, '%s', 'delimiter', '|');
if length(descs) > 0
	fprintf(1, 'OTHER DESCR.:\t');
end
for k = 1:length(descs)
	fprintf(1, '%s', descs{k});
	if k ~= length(descs)
		fprintf(1, '\n\t\t');
	else
		fprintf(1, '\n');
	end
end

tokens = regexp(summary, '<Item Name="Orgname".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, 'ORGANISM:\t%s\n', tokens{1});

tokens = regexp(summary, '<Item Name="Chromosome".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, 'LOCATION:\tChromosome %s, ', tokens{1});

tokens = regexp(summary, '<Item Name="MapLocation".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, 'locus %s, coordinates ', tokens{1});

tokens = regexp(summary, '<Item Name="ChrStart".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, '[%s..', tokens{1});

tokens = regexp(summary, '<Item Name="ChrStop".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
fprintf(1, '%s]\n', tokens{1});

tokens = regexp(summary, '<Item Name="Summary".+?>(.*?)</Item>', 'tokens');
tokens = tokens{1};
if length(tokens{1}) > 0
	fprintf(1, 'SUMMARY:\n%s\n', tokens{1});
else
	fprintf(1, 'SUMMARY:\tN/A\n');
end

return;







function escaped = escape_string(str)
escaped = strrep(str, ' ', '+');
return;

