function methylation = read_illumina_methylation27_samples( ...
	sample_file, probeset_file)

global organism;
genes = organism.Genes;

probesets = read_probesets_illumina_methylation27(probeset_file);

[data, headers] = readtable(sample_file, ...
	'Numeric', {'methylated', 'detection'});

% Map table rows to gene IDs.
valid = probesets.id_to_idx.isKey(data{1});
if any(~valid), error 'Not all probeset IDs were recognized.'; end
ps_idx = cell2mat(probesets.id_to_idx.values(data{1}));

entrez = probesets.entrez(ps_idx);

if rx(headers{2}, 'unmethylated') && rx(headers{3}, 'methylated') && ...
	rx(headers{4}, 'detection')
	
	fprintf('Detected format without beta-values...\n');
	S = (length(headers)-1) / 3;
	for s = 1:S
		tokens = regexp(headers{s*3-1}, '^(.+?)\.Unmethylated', 'tokens');
		if length(tokens) ~= 1, error 'Could not find sample ID.'; end
		token = tokens{1}; methylation.meta.sample_id{s} = token{1};
		
		unmeth(:, s) = data{s*3-1};
		meth(:, s) = data{s*3};
		p_detect(:, s) = data{s*3+1};
	end
	
	methylation.entrez = probesets.entrez(ps_idx);
	methylation.chromosome = probesets.chromosome(ps_idx);
	methylation.position = probesets.position(ps_idx);
	
	methylation.beta = max(unmeth, 0) ./ (max(meth, 0) + max(unmeth, 0) + 100);
else
	error 'Unrecognized file format.';
end

	
	
	
	
	
	
	
	
	
	
	
function probesets = read_probesets_illumina_methylation27(filepath)

global organism;
genes = organism.Genes;

fid = fopen(filepath);

probesets = struct;
probesets.id = {};
probesets.entrez = {};
probesets.chromosome = {};
probesets.position = {};

P = 0;

while 1
	line = fgetl(fid);
	if ~ischar(line), break; end
	
	tokens = textscan(line, '%s', -1, 'Delimiter', ',');
	token = tokens{1};
	if length(token) < 22, continue, end
	
	P = P + 1;
	probesets.id{P,1} = token{1};
	probesets.entrez{P,1} = token{22}(8:end);
	probesets.chromosome{P,1} = token{9};
	probesets.position{P,1} = token{10};
end

fclose(fid);

probesets.entrez = str2double(probesets.entrez);
valid = ~strcmp(probesets.chromosome, 'Chr');
probesets = filter_struct(probesets, valid);
probesets.chromosome = chromosome_sym2num(probesets.chromosome);
probesets.position = str2double(probesets.position);

probesets = liftover(probesets, 'hg18->hg19');

probesets.id_to_idx = containers.Map(probesets.id, ...
	num2cell(1:length(probesets.id)));

