
function [] = vcf_discard_contingent(vcf_file, alpha)

fprintf('Reading VCF file...\n');
[data, headers] = readtable(vcf_file);

V = length(data{1});

for k = 1:length(headers)
	if ~iscellstr(data{k}), continue, end
	if any(rx(data{k}(1:min(V, 100)), '[01]/[01]:'))
		first_sample_col = k; break;
	end
end

S = length(headers) - first_sample_col + 1;

samples = headers(first_sample_col:end)';
alt_reads = nan(V, S);
total_reads = nan(V, S);

for s = 1:S
	info = data{first_sample_col + s - 1};
	str = sprintf('%s\n', info{:});
	sdata = textscan(str, '%*s %d %d %*[^\n]', 'Delimiter', ':,', ...
		'ReturnOnError', false);
	alt_reads(:, s) = sdata{1};
	total_reads(:, s) = sdata{2};
end

data = cat(2, data{:});

out = fopen('contingent.vcf', 'w');
fprintf(out, '%s\t', headers{1:end-1});
fprintf(out, '%s\n', headers{end});

for v = 1:V
	if mod(v, 100) == 0, fprintf('%d / %d\n', v, V); end

	nz = (total_reads(v, :) ~= 0);
	p = chi2cont([alt_reads(v,nz); total_reads(v,nz) - alt_reads(v,nz)]);
	if p > alpha, continue, end
	
	fprintf(out, '%s\t', data{v, 1:end-1});
	fprintf(out, '%s', data{v, end});
	if alpha == 1, fprintf(out, '\t%.2f', -log10(p)); end
	fprintf(out, '\n');
end

fclose(out);
