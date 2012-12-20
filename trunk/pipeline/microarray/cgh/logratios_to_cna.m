
function cna = logratios_to_cna(logratios, gender)

global organism;

if ischar(gender)
	gender = repmat({gender}, 1, length(logratios.meta.sample_id));
end

gender_num = nan(1, length(gender));
gender_num(rx(gender, '^male')) = 1;
gender_num(rx(gender, 'female')) = 2;

ploidy = organism.Chromosomes.Ploidy(:, gender_num);

if isfield(logratios, 'rows') && isfield(logratios.rows, 'chromosome')
	
	chr = logratios.rows.chromosome;
	valid = ~isnan(chr);
	if any(~valid)
		logratios = filter_rows(logratios, valid);
		fprintf('Discarded %d features with unknown chromosome.\n',sum(~valid));
	end

	cna = logratios;
	cna.mean = ploidy(chr, :) .* 2.^logratios.mean - ploidy(chr, :);
	cna.mean = max(cna.mean, -ploidy(chr, :));
	cna.mean(isnan(logratios.mean)) = NaN;   % max() removes NaNs

elseif rx(logratios.meta.type, 'Gene copy.*log')
	
	[P, S] = size(logratios.mean);
	
	gidx = gene_idx(logratios.rows.gene_symbol);
	valid = ~isnan(gidx);
	if any(~valid)
		logratios = filter_rows(logratios, valid);
		gidx = gidx(valid);
		fprintf('Discarded %d genes with an unknown name.\n', sum(~valid));
	end
	
	chr = organism.Genes.Chromosome(gidx);
	valid = ~isnan(chr);
	if any(~valid)
		logratios = filter_rows(logratios, valid);
		gidx = gidx(valid);
		fprintf('Discarded %d genes due to unknown chromosome.\n', sum(~valid));
	end
	
	cna = logratios;
	cna.mean = ploidy(chr, :) .* 2.^logratios.mean - ploidy(chr, :);
	cna.mean = max(cna.mean, -ploidy(chr, :));
	cna.mean(isnan(logratios.mean)) = NaN;
	
else
	error 'Cannot recognize logratio data format.';
end

