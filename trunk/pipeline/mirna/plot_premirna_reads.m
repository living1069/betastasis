
% Author: Matti Annala <matti.annala@tut.fi>

function [] = plot_premirna_reads(reads, varargin)
	
global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

hide_mature_reads = false;

S = length(reads.Raw);

for s = 1:S
	fprintf(1, 'Plotting reads for sample %s...\n', ...
		reads.Meta.Sample.ID{s});
	
	trim_len = 18;
	trimmed = trim_reads(filter_query(reads, s), trim_len);
	
	al = align_reads(trimmed, 'pre_mirnas', 'MaxMismatches', 2, ...
		'Columns', 'target,offset', varargin{:});
			
	for p = 1:length(pre_mirnas.Name)
		pre_reads = strcmp(pre_mirnas.Name{p}, al.Target);
		
		offsets = al.Offset(pre_reads);
		
		if hide_mature_reads
			mature_reads = false(size(offsets));
			for m = 1:pre_mirnas.MatureCount(p)
				idx = pre_mirnas.Matures(p, m);
				mature_reads(offsets >= pre_mirnas.MatureOffsets(p, m) & ...
					offsets + trim_len <= pre_mirnas.MatureOffsets(p, m) + ...
					length(mirnas.Sequence{idx})) = true;
			end
			offsets = offsets(~mature_reads);
		end
		
		read_count_dist = zeros(1, length(pre_mirnas.Sequence{p}) + 1);
		[length(read_count_dist), max(offsets) + trim_len-1]
		for m = 1:length(offsets)
			span = offsets(m):offsets(m)+trim_len-1;
			read_count_dist(span) = read_count_dist(span) + 1;
		end
		
		figure; bar(1:length(read_count_dist), read_count_dist);
		title(sprintf('%s in sample %s', pre_mirnas.Name{p}, ...
			reads.Meta.Sample.ID{s}), 'Interpreter', 'none');
		xlabel('pre-miRNA sequence offset');
		ylabel('Number of overlapping reads');
		saveas(gcf, ['~/pre_mirnas/' pre_mirnas.Name{p} '.pdf']);
	end
end

