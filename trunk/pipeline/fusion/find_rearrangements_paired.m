
% Author: Matti Annala <matti.annala@tut.fi>

function rearrangements = find_rearrangements_paired(reads, varargin)

global organism;

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	%if strcmpi(varargin{k}, 'PriorRearrangements')
	%	prior_rearrangements = varargin{k+1};
	%	drop_args(k:k+1) = true;
	%	continue;
	%end
end
varargin = varargin(~drop_args);

S = length(reads.Raw);

rearrangements.Fusions = cell(1, S);
for s = 1:S
	rearrangements.Fusions{s} = find_fusions(filter_query(reads, s));
end

rearrangements.Meta = reads.Meta;
rearrangements.Meta.Type = 'Genetic rearrangements';
rearrangements.Meta.Organism = organism.Name;
	









function fusions = find_fusions(reads)
	
global organism;
transcripts = organism.Transcripts;
exons = organism.Exons;

al = align_reads(reads, 'transcripts', 'IgnorePairs', true, ...
	'AllowAlignments', 20, 'ReportAlignments', 20, 'MaxMismatches', 1, ...
	'Columns', 'read,target');

read_ids = al.ReadID;
tx_indices = transcript_idx(al.Target);




fprintf(1, 'Searching for rearrangement events...\n');
run_ends = [ find(strcmp(read_ids(1:end-1), read_ids(2:end)) == 0); ...
	length(read_ids) ];
run_lengths = diff([0; run_ends]);
run_starts = run_ends - run_lengths + 1;

read_id_to_run = containers.Map(read_ids(run_starts), ...
	num2cell(1:length(run_lengths)));

found_fusions = 0;
fusion_map = containers.Map;

fusions = struct;
fusions.Exons = zeros(0, 2);
fusions.ReadCount = [];

progress = Progress;

for r = 1:length(run_starts)
	id = read_ids{run_starts(r)};
	
	if ~strcmp('/2', id(end-1:end)), continue, end
	
	if read_id_to_run.isKey([id(1:end-2) '/1'])
		l = read_id_to_run([id(1:end-2) '/1']);
		left_run = run_starts(l):run_ends(l);
		right_run = run_starts(r):run_ends(r);
		
		left_txs = tx_indices(left_run);
		right_txs = tx_indices(right_run);
		
		left_genes = transcripts.Gene(left_txs);
		right_genes = transcripts.Gene(right_txs);
		
		% Our goal is to try to find the simplest hypothesis that can explain
		% the origin of each read that did not align to the transcriptome.
		if length(left_genes) > 1 || length(right_genes) > 1
			continue;
			
		elseif ~isempty(intersect(left_genes, right_genes))
			
			% We have an alignment pair that indicates that the tags come from
			% the same gene. Hence we need not hypothesize a fusion event.
			continue;
		else
			%left_strands = al.Strand(left_run);
			%right_strands = al.Strand(right_run);
			
			%left_sequences = al.Sequence(left_run);
			%right_sequences = al.Sequence(right_run);
			
			sorted_genes = sort([left_genes(1) right_genes(1)]);
			key = sprintf('%d,%d', sorted_genes(1), sorted_genes(2));
			
			if ~fusion_map.isKey(key)
				found_fusions = found_fusions + 1;
				fusion_map(key) = found_fusions;
				fusions.Genes(found_fusions, :) = sscanf(key, '%d,%d');
				fusions.ReadCount(found_fusions, 1) = 1;
				%fusions.ReadSequences{found_fusions} = ...
				%	{ left_sequences{k}, right_sequences{m} };
			else
				idx = fusion_map(key);
				fusions.ReadCount(idx) = fusions.ReadCount(idx) + 1;
				%fusions.ReadSequences{idx}(end+1, :) = ...
				%	{ left_sequences{k}, right_sequences{m} };
			end
		end
	end
	
	progress.update(r / length(run_starts));
end

fprintf(1, 'Found %d fusion events.\n', found_fusions);


