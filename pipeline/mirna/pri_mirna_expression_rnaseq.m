
% Author: Matti Annala <matti.annala@tut.fi>

function expr = pri_mirna_expression_rnaseq(reads, varargin)

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

pri_mirnas = find_pri_mirna();

normalization = 'none';

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Normalization')
		if strcmpi('none', varargin{k+1})
			normalization = 'none';
		elseif strcmpi('rpkm', varargin{k+1})
			normalization = 'rpkm';
		else
			error(['Invalid normalization method specified. Valid ' ...
				   'alternatives include "none" and "rpkm".']);
		end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = length(reads.Raw);

pri_mirnas.Name = organism.pre_miRNA.Name;
primir_name_to_idx = containers.Map(pri_mirnas.Name, ...
	num2cell(1:length(pri_mirnas.Name)));
	
mirna_lengths = zeros(length(mirnas.Sequence), 1);
for k = 1:length(mirna_lengths)
	mirna_lengths(k) = length(mirnas.Sequence{k});
end

expr = struct;
expr.Mean = zeros(length(pri_mirnas.Name), S);

expr.Meta = struct;
if isstruct(reads), expr.Meta = reads.Meta; end

expr.Meta.Type = 'pri-miRNA expression';
expr.Meta.Organism = organism.Name;
expr.Meta.DateCreated = repmat({datestr(now)}, S, 1);
expr.Meta.Normalization = repmat({normalization}, S, 1);
expr.Meta.TotalSeqReads = zeros(S, 1);

for s = 1:length(reads.Raw)
	extracted = extract_reads(filter_query(reads, s));
	
	fprintf(1, 'Aligning reads to miRBase using Bowtie...\n');
	al = align_reads(extracted, pri_mirnas, ...
		'MaxMismatches', 2, 'AllowAlignments', 10, ...
		'Columns', 'target,offset', varargin{:});
	
	expr.Meta.TotalSeqReads(s) = al.TotalReads;

	targets = cell2mat(primir_name_to_idx.values(al.Target));
	read_offsets = al.Offset;

	fprintf(1, 'Summarizing pri-miRNA expression levels...\n');

	for t = targets'
		expr.Mean(t, s) = expr.Mean(t, s) + 1;
	end
	
	if strcmp('none', normalization)
		continue;
	elseif strcmp(normalization, 'rpkm')
		fprintf(1, 'Normalizing pri-miRNA expression levels using RPM...\n');
		expr.Mean(:, s) = expr.Mean(:, s) / (al.TotalReads / 1e6);
	end
end

