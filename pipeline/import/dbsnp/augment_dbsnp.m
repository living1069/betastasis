function [] = augment_dbsnp(vcfa_file)

[status, out] = unix(sprintf( ...
	'~/tools/bedtools/bin/intersectBed -wao -a <(cut ) -b <()', ...
	vcfa_file));

%[status, ~] = unix( ...
%	'~/pipeline/import/dbsnp/dbsnp_flat_to_tabular.py');
%if status ~= 0, error 'dbSNP import failed.'; end

fifo = 'dbsnp.txt.fifo';
unix(['mkfifo ' fifo]);
unix(sprintf('gunzip -c dbsnp*.gz > %s &', fifo));
[data, headers] = readtable(fifo);

dbsnp.chromosome = chromsome_sym2num(data{1});
dbsnp.position = str2double(data{2});
dbsnp.id = data{3};
dbsnp.alleles = data{4};
dbsnp.heterozygosity = data{5};

bad = isnan(dbsnp.position);
dbsnp.chromosome(bad)
dbsnp.position(bad)
	
%dbsnp.map = containers.Map(strcat(data{1}, ':', data{2}), ...
%	num2cell(1:length(data{1})));

