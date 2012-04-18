function [] = build_bowtie2_indices()

global organism;

index_root = [ppath '/tools/bowtie2/indexes/' flatten_str(organism.Name) ...
	'/' flatten_str(organism.Version)];
[~, ~] = mkdir(index_root);

bowtie_build = [ppath '/tools/bowtie2/bowtie2-build'];
fasta_tmp = ptemp;


write_seq_fasta(organism.Transcripts, fasta_tmp);
[status, out] = unix(sprintf('%s %s %s/transcripts', ...
	bowtie_build, fasta_tmp, index_root));
if status ~= 0, fprintf(1, 'bowtie2_build failed:\n%s\n', out); end

write_seq_fasta(organism.Exons, fasta_tmp);
status = unix(sprintf('%s %s %s/exons', ...
	bowtie_build, fasta_tmp, index_root));
if status ~= 0, fprintf(1, 'bowtie2_build failed:\n%s\n', out); end


