function [color, quality] = seq_read_type(reads)

file = fopen(reads);

quality = 0;
color = -1;

while 1
	line = fgetl(file);
	if line == -1, break, end
	
	if line(1) == '#' || line(1) == '>', continue, end
	
	if line(1) == '@', quality = 1; end

	if regexpi(line, '^[tcgarykmswbdhvn]+$')
		color = 0; break;
	elseif regexpi(line, '^[ACTG][0123\.]+$')
		color = 1; break;
	else
		error 'Read sequences are in an unrecognized format.';
	end
end

fclose(file);

