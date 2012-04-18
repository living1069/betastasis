function [] = sms_to_fasta_batch(sms_reads, fasta_dataset)

global pipeline_config;

setenv('HELICOS_TMPDIR', pipeline_config.TempDir);
pwd = cd;
cd(pipeline_config.TempDir);

fprintf(1, 'Converting SMS reads to FASTA...\n');
progress = Progress;

seq_files = seq_resource_files(sms_reads);

[~, ~] = mkdir(raw_dir);

for k = 1:length(sms_reads.SequenceResource)
	fasta_reads = [raw_dir '/' regexprep(sms_reads.SequenceResource{k}, ...
		'^.*/([^/]*).sms$', '$1', 'ignorecase')];
	status = unix(sprintf(['%s/tools/helisphere/bin/sms2txt --fasta ' ...
		'--input_file %s --output_file_prefix %s'], ...
		ppath, seq_files{k}, fasta_reads));
	if status ~= 0, error 'SMS to FASTA conversion failed.'; end
	
	progress.update(k / length(sms_reads.SequenceResource));
end

cd(pwd);

