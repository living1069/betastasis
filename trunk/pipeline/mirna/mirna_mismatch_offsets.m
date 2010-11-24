function counts = mirna_mismatch_offsets(alignments)

fid = fopen(alignments);
data = textscan(fid, '%*s %*s %*s %d %*s %s %*d');
fclose(fid);

offsets = data{1};
qualities = data{2};

perfect_count = 0;
total_count = length(qualities);
counts = zeros(1, 100);

for k = 1:length(qualities)
	bad = find(qualities{k} == '!');
	if isempty(bad)
		perfect_count = perfect_count + 1;
		continue;
	end
	
	counts(offsets(k) + bad(1)) = counts(offsets(k) + bad(1)) + 1;
end

perfect_count
total_count
