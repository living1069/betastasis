
% Author: Matti Annala <matti.annala@tut.fi>

function probesets = read_probesets_illumina_methylation27(filepath)

global organism;
genes = organism.Genes;

fid = fopen(filepath);

probesets = struct;
probesets.id = cell(1e5, 1);
probesets.entrez = nan(1e5, 1);

P = 0;

while 1
	line = fgetl(fid);
	if ~ischar(line), break; end
	
	tokens = regexp(line, '^(cg.+?),.*,GeneID:(\d+),', 'tokens');
	if length(tokens) ~= 1, continue, end
	
	token = tokens{1}; P = P + 1;
	probesets.id{P} = token{1};
	probesets.entrez(P) = str2double(token{2});
end

probesets.id = probesets.id(1:P);
probesets.entrez = probesets.entrez(1:P);

