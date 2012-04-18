function [] = print_chimeras(chimeras, varargin)

global organism;

min_frequency = 0;
min_anchor_len = 20;
min_read_anchor_len = 10;
negative_samples = [];
max_avg_mismatches = Inf;
ranking = 'chisquare';

% Organism specific genes that are blacklisted by default.
if strcmpi(organism.Name, 'Homo sapiens')
	blacklist = { 'GCGCTTGA' };
end

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Blacklist')
		blacklist = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MinFrequency')
		min_frequency = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'max.*(avg|average).*mismatch')
		max_avg_mismatches = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'MinAnchorLen')
		min_anchor_len = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'MinReadAnchorLen')
		min_read_anchor_len = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'NegativeSamples')
		negative_samples = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

C = length(chimeras.sequence);
S = size(chimeras.reads, 2);

fprintf(1, '%d total reads aligned to %d chimeric junctions...\n', ...
	sum(sum(cellfun(@length, chimeras.reads))), C);




fprintf(1, 'Filtering chimeras according to specified criteria...\n');

discard = zeros(1, C);
for c = 1:C
	
	% Filter out chimeras that start with RNY5.
	if ~isempty(blacklist)
		reads = cat(1, chimeras.reads{c, :});
		blacklisted_reads = 0;
		for r = 1:length(reads)
			for b = 1:length(blacklist)
				if strfind(reads{r}, blacklist{1}) == 1
					discard(c) = 3; continue;
				end
			end
		end
	end
	
	
	% Check that the reads do not have an aberrantly high number of mismatches.
	if max_avg_mismatches < Inf
		mismatches = [0 0];
		reads = cat(1, chimeras.reads{c, :});
		for r = 1:length(reads)
			pos = strfind(chimeras.sequence{c}, reads{r});
			mismatches = mismatches + [sum(seq ~= ref_seq), 1];
		end
		
		if mismatches(1) / mismatches(2) > max_avg_mismatches
			discard(c) = 1; continue;
		end
	end
	
	if ~isempty(negative_samples)
		if any(~cellfun(@isempty, chimeras.reads(c, negative_samples)))
			discard(c) = 2; continue;
		end
	end
end


fprintf(1, '%d / %d (%.1f%%) discarded due to too many mismatches.\n', ...
	sum(discard == 1), length(discard), ...
	100 * sum(discard == 1) / length(discard));
fprintf(1, '%d / %d (%.1f%%) discarded due to presence in controls.\n', ...
	sum(discard == 2), length(discard), ...
	100 * sum(discard == 2) / length(discard));


fprintf(1, '%d / %d (%.1f%%) chimeras passed all criteria.\n', ...
	sum(discard == 0), C, sum(discard == 0) / C * 100);

	
chimeras.sequence = chimeras.sequence(discard == 0);
chimeras.reads = chimeras.reads(discard == 0, :);
C = length(chimeras.sequence);


% Next we sort the chimeras according to their chi-square deviation from
% the uniform distribution.
if rx(ranking, 'chi.*sq')
	chisq = nan(1, C);
	for c = 1:C
		sr = cellfun(@length, chimeras.reads(c, :));
		expected_reads = sum(sr) / S;
		chisq(c) = sum((sr - ones(1, S) * expected_reads).^2 / expected_reads);
	end
	
	[~, order] = sort(chisq, 'descend');
end



% Finally we print the remaining candidates.
fid = fopen('~/chimeras.txt', 'W');
for c = order
	fprintf(fid, 'Chimera %s:\n', chimeras.sequence{c});
	for s = 1:S
		for r = 1:length(chimeras.reads{c, s})
			fprintf(fid, '- %d: %s\n', s, chimeras.reads{c, s}{r});
		end
	end
end
fclose(fid);
