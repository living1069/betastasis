function fasta_reads = sms_to_fasta(sms_reads, varargin)

global pipeline_config;

for k = 1:2:length(varargin)
	error('Unrecognized option "%s".', varargin{k});
end

setenv('HELICOS_TMPDIR', pipeline_config.TempDir);

fprintf(1, 'Filtering SMS reads...\n');
filtered_sms = ptemp;
status = unix(sprintf(['%s/tools/helisphere/bin/filterSMS ' ...
	'--input_file %s --output_file %s --minlen 24 --quality 10 ' ...
	'--minscore 4'], ppath, sms_reads, filtered_sms));
	
fprintf(1, 'Converting SMS reads to FASTA...\n');
fasta_reads = ptemp;
[status, ~] = unix(sprintf(['%s/tools/helisphere/bin/sms2txt --fasta ' ...
	'--input_file %s --output_file_prefix %s'], ...
	ppath, filtered_sms, fasta_reads));
if status ~= 0, error 'SMS to FASTA conversion failed.'; end
	
safe_delete(filtered_sms);
fasta_reads = [fasta_reads '.fa'];

