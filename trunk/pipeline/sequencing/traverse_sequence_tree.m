function [] = traverse_sequence_tree(reads, seed, max_depth)

paths = { seed };
read_len = 34;

fasta_tmp = ptemp();
index_tmp = ptemp();
alignments_tmp = ptemp();

while 1
	if isempty(paths), break, end
	path = paths{end};
	paths = paths(1:end-1);
	
	if length(path) > max_depth, continue, end
	
	fid = fopen(fasta_tmp, 'w');
	for nuc = 'TCGA'
		fprintf(fid, '>%s\n%s%s\n', nuc, path, nuc);
	end
	fclose(fid);
	
	[status, out] = unix([ppath '/tools/bowtie/bowtie-build -C ' ...
		fasta_tmp ' ' index_tmp]);
	if status ~= 0, error('%s\n', out); end

	trim = abs(length(path) - read_len);
	[status, out] = unix(sprintf('%s/tools/bowtie/bowtie -C -f -p4 -v0 -k4 --trim3 %d --suppress 5,6,7,8 %s %s > %s', ppath, trim, index_tmp, reads, alignments_tmp));

	%fprintf(1, '%s\n', out);
	if status ~= 0, error('%s\n', out); end
	
	fid = fopen(alignments_tmp);
	data = textscan(fid, '%*s %*s %s %*s');
	fclose(fid);
	
	read_nucs = data{1};
	for nuc = 'TCGA'
		ratio = sum(strcmp(nuc, read_nucs)) / length(read_nucs);
		%fprintf(1, '%s: %.1f\n', nuc, ratio * 100);
		if ratio > 0.30
			fprintf(1, '%s (%.1f%%)\n', [path nuc], ratio * 100);
			paths{end + 1} = [path nuc];
		end
	end
end


