function chr_lengths = read_chromosomes(filepath)

chr_lengths = zeros(25, 1);

chromosomes = {};
for k = 1:22, chromosomes{k} = ['chr' num2str(k)]; end
chromosomes{23} = 'chrX';
chromosomes{24} = 'chrY';
chromosomes{25} = 'chrM';

fid = fopen(filepath);
if fid == -1
	fprintf(1, 'File %s could not be opened.', filepath);
	return;
end

chr = '';
seq = repmat(' ', 1, 3e8);

while 1
	line = fgetl(fid);
	if line == -1
		chrnum = find(strcmp(chr, chromosomes));
		sequence = seq(1:pos-1);
		chr_lengths(chrnum) = length(sequence);
		save(chr, 'sequence');
		break;
	end
		
	if line(1) == '>'
		if ~isempty(chr)
			chrnum = find(strcmp(chr, chromosomes));
			sequence = seq(1:pos-1);
			chr_lengths(chrnum) = length(sequence);
			save(chr, 'sequence');
		end
		chr = line(2:end);
		pos = 1;
		fprintf(1, 'Reading the sequence for chromosome %s...\n', chr);
	end
	
	len = length(line);
	seq(pos:pos+len-1) = line;
	pos = pos + len;
end

