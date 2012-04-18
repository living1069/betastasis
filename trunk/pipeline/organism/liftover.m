function features = liftover(features, direction)

global organism;
chromosomes = organism.Chromosomes;

old_file = ptemp;
new_file = ptemp;
unmapped_file = ptemp
fid = fopen(old_file, 'W');
for k = 1:length(features.chromosome)
	fprintf(fid, '%s\t%d\t%d\n', ...
		['chr' chromosomes.Name{features.chromosome(k)}], ...
		features.position(k) - 1, features.position(k));
end
fclose(fid);

chain_file = '';
if regexpi(direction, 'hg18\s*->\s*hg19')
	chain_file = sprintf('%s/tools/liftover/hg18ToHg19.over.chain', ppath);
elseif regexpi(direction, 'hg19\s*->\s*hg18')
	chain_file = sprintf('%s/tools/liftover/hg19ToHg18.over.chain', ppath);
else
	error 'Unknown chain.';
end

status = unix(sprintf('%s/tools/liftover/liftOver %s %s %s %s', ...
	ppath, old_file, chain_file, new_file, unmapped_file));
if status ~= 0, error 'Liftover failed.'; end

fid = fopen(new_file);
data = textscan(fid, '%s %d %d', 'Delimiter', '\t');
fclose(fid);

fid = fopen(unmapped_file);
undata = textscan(fid, '%s %d %d', 'Delimiter', '\t', 'CommentStyle', '#');
fclose(fid);

if length(data{1}) + length(undata{1}) ~= length(features.chromosome)
	error 'Mapped + unmapped count does not match input count.';
end

undata{1} = chromosome_sym2num(undata{1});

unmapped = false(length(features.chromosome), 1);
for u = 1:length(undata{1})
	un = find(features.chromosome == undata{1}(u) & ...
		features.position == undata{3}(u));
	if length(un) ~= 1
		error 'Cannot find unique match for unmapped entry.';
	end
		
	unmapped(un) = true;
end

features.chromosome(~unmapped) = chromosome_sym2num(data{1})';
features.position(~unmapped) = data{3};

features.chromosome(unmapped) = NaN;
features.position(unmapped) = NaN;

safe_delete(old_file);
safe_delete(new_file);
safe_delete(unmapped_file);

