function [] = darned_to_annovar_generic(darned_file)

[data, headers] = readtable(darned_file);
chr = data{strcmpi(headers, 'chr')};
offset = str2double(data{strcmpi(headers, 'coordinate')});
strand = data{strcmpi(headers, 'strand')};
ref_allele = data{strcmpi(headers, 'inchr')};
obs_allele = data{strcmpi(headers, 'inrna')};

obs_allele = strrep(obs_allele, 'I', 'G');

% Complement the alleles if we are talking about the opposite strand.
minus = strcmp(strand, '-');
for k = find(minus)'
	ref_allele{k} = seqcomplement(ref_allele{k});
	obs_allele{k} = seqcomplement(obs_allele{k});
end

bed_file = regexprep(darned_file, '\.txt$', '.annovar.txt');

fid = fopen(bed_file, 'W');
for k = 1:length(chr)
	fprintf(fid, '%s\t%d\t%d\t%s\t%s\n', chr{k}, offset(k), offset(k), ...
		ref_allele{k}, obs_allele{k});
end
fclose(fid);

