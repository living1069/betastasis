function [] = entrez_transcript(transcript_names, organism)

organism_name = 'Homo sapiens';
if nargin == 2
	organism_name = organism;
end

if iscellstr(transcript_names)
	for k = 1:length(transcript_names)
		summarize_transcript(transcript_names{k}, organism_name);
		if k ~= length(transcript_names)
			fprintf(1, ['==============================================' ...	
			            '===============================\n']);
			pause(3);
		end
	end
elseif ischar(transcript_names)
	summarize_transcript(transcript_names, organism_name);
else
	error 'Transcript names parameter has an invalid type.';
end

return;




function [] = summarize_transcript(transcript_name, organism_name)

entrez_base_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils';
database = 'nuccore';

search_results = urlread(sprintf('%s/esearch.fcgi?db=%s&term=%s+AND+(%s[Organism])', ...
	entrez_base_url, database, transcript_name, escape_string(organism_name)));

ids = regexp(search_results, '<Id>(\w+)</Id>', 'tokens');
if length(ids) == 0
	fprintf(1, 'No matches found.\n');
	return;
end

best_match = ids{1}; best_match = best_match{1};

info = urlread(sprintf('%s/efetch.fcgi?db=%s&retmode=xml&id=%s', ...
	entrez_base_url, database, best_match));

%info_query = sprintf('%s/einfo.fcgi?db=%s&id=%s', entrez_base_url, ...
%	database, best_match);

%info = urlread(info_query);

tokens = regexp(info, '<GBSeq_accession-version>(.*?)</', 'tokens');
tokens = tokens{1};
fprintf(1, 'ACCESSION:\t%s\n', tokens{1});

tokens = regexp(info, '<GBSeq_definition>(.*?)</', 'tokens');
tokens = tokens{1};
fprintf(1, 'DESCRIPTION:\t%s\n', tokens{1});

tokens = regexp(info, '<GBSeq_organism>(.*?)</', 'tokens');
tokens = tokens{1};
fprintf(1, 'ORGANISM:\t%s\n', tokens{1});

tokens = regexp(info, '<GBSeq_length>(.*?)</', 'tokens');
tokens = tokens{1};
fprintf(1, 'LENGTH:\t\t%s bp\n', tokens{1});

tokens = regexp(info, '<GBSeq_create-date>(.*?)</', 'tokens');
tokens = tokens{1};
fprintf(1, 'CREATED:\t%s\n', tokens{1});

tokens = regexp(info, '<GBSeq_update-date>(.*?)</', 'tokens');
tokens = tokens{1};
fprintf(1, 'UPDATED:\t%s\n', tokens{1});

tokens = regexp(info, '<GBSeq_comment>(.*?)</', 'tokens');
tokens = tokens{1};
comment = tokens{1};
comment = strrep(comment, '&apos', '''');
fprintf(1, '\n%s\n', comment);


tokens = regexp(info, '<GBReference_title>(.*?)</', 'tokens');
if length(tokens) > 0
	fprintf(1, '\nREFERENCE TITLES:\n');
end
for k = 1:length(tokens)
	token = tokens{k};
	fprintf(1, '- %s\n', token{1});
end

return;







function escaped = escape_string(str)
escaped = strrep(str, ' ', '+');
return;

