function fasta = read_fasta(fasta_file)

fid = fopen(fasta_file);

fasta = struct;
fasta.name = {};
fasta.sequence = {};

k = 0;

while 1
	line = fgetl(fid);
	if ~ischar(line), break, end
		
	if line(1) == '#', continue, end
	
	if line(1) == '>'
		k = k + 1;
		fasta.name{k, 1} = line(2:end);
		fasta.sequence{k, 1} = '';
	else
		fasta.sequence{k} = [fasta.sequence{k}, line];
	end
end

fclose(fid);

