function refgene = read_ucsc_refgene(filepath)

global organism;
genes = organism.Genes;

refgene = struct;
refgene.Transcripts = struct;
refgene.Transcripts.Name = {};
refgene.Transcripts.Chromosome = [];
refgene.Transcripts.Exons = {};
refgene.Transcripts.CDS = [];

fid = fopen(filepath);
data = textscan(fid, '%s %s %s %s %d %d %d %d %d %s %s %s %s %s %s %s', ...
	'Delimiter', '\t');
fclose(fid);

transcript = data{2};
chromosome = data{3};
strand = data{4};
left = data{5};
right = data{6};
cds_left = data{7};
cds_right = data{8};
exon_starts = data{10};
exon_ends = data{11};

drop = false(length(transcript), 1);
for ts = 1:length(transcript)
	refgene.Transcripts.Name{end+1, 1} = transcript{ts};
	refgene.Transcripts.Chromosome(end+1, 1) = chromosome_sym2num( ...
		chromosome{ts});
	
	estart = textscan(exon_starts{ts}, '%d,', -1);
	eend = textscan(exon_ends{ts}, '%d,', -1);
	refgene.Transcripts.Exons{end+1, 1} = [estart{1} eend{1}];
	
	refgene.Transcripts.CDS(end+1, :) = [cds_left(ts) cds_right(ts)];
	
	if drop(ts), continue, end
	
	dups = find(strcmpi(transcript{ts}, transcript));
	if length(dups) > 1
		drop(dups) = true;

		% If all but one transcript have been annotated with a strange
		% chromosome, we discard all the anomalous transcripts. If there are
		% multiple non-anomalous transcripts, then we keep them all.
		weird_chr = false(length(dups), 1);
		for k = 1:length(dups)
			if regexpi(chromosome{dups(k)}, '^chr\d+_.*$|chrUn')
				weird_chr(k) = true;
			end
		end
		dups = dups(~weird_chr);
		
		if length(dups) == 1
			drop(dups(1)) = false;
			continue;
		elseif length(dups) == 0
			continue;
		end
		
		% Enable this code if you want to see a trace of duplicated transcripts.
		if 0
			fprintf(1, 'Transcript %s, %d dups:\n', transcript{ts}, ...
				length(dups));
			
			for k = 1:length(dups)
				fprintf(1, '- %s (%s), [%d..%d]\n', chromosome{dups(k)}, ...
					strand{dups(k)}, left(dups(k)), right(dups(k)));
			end
			
			for k = 2:length(dups)
				continue;
				if left(dups(1)) ~= left(dups(k))
					error 'Left coordinate mismatch.';
				end
				if right(dups(1)) ~= right(dups(k))
					error 'Right coordinate mismatch.';
				end
				if ~strcmp(strand{dups(1)}, strand{dups(k)})
					error 'Strand mismatch.';
				end
				if ~strcmp(chromosome{dups(1)}, chromosome{dups(k)})
					error 'Chromosome mismatch.';
				end
			end
		end
	end
end

transcript = transcript(~drop);
chromosome = chromosome(~drop);
strand = strand(~drop);
left = left(~drop);
right = right(~drop);

transcript_to_idx = containers.Map(transcript, num2cell(1:length(transcript)));

refgene.Chromosome = nan(length(genes.Name), 1);
refgene.Position = nan(length(genes.Name), 2);
refgene.Strand = repmat(' ', length(genes.Name), 1);

for g = 1:length(genes.Name)
	gene_left = Inf;
	gene_right = -Inf;
	
	gene_ts = organism.Transcripts.Name(genes.Transcripts( ...
		g, 1:genes.TranscriptCount(g)));
	found = transcript_to_idx.isKey(gene_ts);
	
	ts_indices = cell2mat(transcript_to_idx.values(gene_ts(found)));
	for k = 1:length(ts_indices)
		idx = ts_indices(k);
		if left(idx) < gene_left, gene_left = left(idx); end
		if right(idx) > gene_right, gene_right = right(idx); end
	end
	
	chr = unique(chromosome(ts_indices));
	if length(chr) ~= 1, continue, end
	
	str = unique(strand(ts_indices));
	if length(str) ~= 1, continue, end
	
	refgene.Chromosome(g) = chromosome_sym2num(chr{1});
	refgene.Strand(g, 1) = str{1};
	refgene.Position(g, :) = [gene_left gene_right];
end

