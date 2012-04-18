function index_path = bowtie_index(feature)

global organism;

feature = lower(feature);
if regexp(feature, '^(genome|transcripts|exons|mirnas|pre_mirnas)$')
	index_path = feature;
else
	error('Bowtie index requested for unsupported feature "%s".', feature);
end

index_path = [ppath '/tools/bowtie/indexes/' flatten_str(organism.Name) ...
	'/' flatten_str(organism.Version) '/' index_path];

