function [] = build_bowtie_indices()

global organism;

index_root = [ppath '/tools/bowtie/indexes/' organism.Version];
bowtie_build = [ppath '/tools/bowtie/bowtie-build'];
fasta_tmp = ptemp;



if isfield(organism, 'Transcripts')
	write_seq_fasta(organism.Transcripts, fasta_tmp);
	status = unix(sprintf('%s %s %s/transcripts', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end

	status = unix(sprintf('%s -C %s %s/transcripts_colorspace', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end
end	

if isfield(organism, 'Exons')
	write_exons_fasta(fasta_tmp);
	status = unix(sprintf('%s %s %s/exons', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end

	status = unix(sprintf('%s -C %s %s/exons_colorspace', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end
end

if isfield(organism, 'miRNA')
	write_seq_fasta(organism.miRNA, fasta_tmp);
	status = unix(sprintf('%s %s %s/mirnas', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end

	status = unix(sprintf('%s -C %s %s/mirnas_colorspace', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end
end

if isfield(organism, 'pre_miRNA')
	write_seq_fasta(organism.pre_miRNA, fasta_tmp);
	status = unix(sprintf('%s %s %s/pre_mirnas', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end

	status = unix(sprintf('%s -C %s %s/pre_mirnas_colorspace', ...
		bowtie_build, fasta_tmp, index_root));
	if status ~= 0, return, end
end

