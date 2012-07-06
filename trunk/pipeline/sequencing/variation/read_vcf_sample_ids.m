
function sample_ids = read_vcf_sample_ids(vcf_file)

[data, headers] = readtable(vcf_file, 'Comment', '^##', 'NumLines', 1);

for k = 1:length(headers)
	if rx(data{k}{1}, '[01]/[01]:')
		first_sample_col = k;
		break;
	end
end

sample_ids = headers(first_sample_col:end)';

