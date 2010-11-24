function [color, quality] = seq_read_type(reads)

file = fopen(reads, 'r');

quality = -1;
color = -1;

while 1
	line = fgetl(file);
	if line == -1, break, end
	
	if strcmp('#', line(1)), continue, end
	
	if quality ~= -1
		if regexpi(line, '^[tcgarykmswbdhvn]+$')
			color = 0;
			break;
		elseif regexpi(line, '^[ACTG][0123\.]+$')
			color = 1;
			break;
		else
			error 'Read sequences are in an unrecognized format.';
		end
	else
		if strcmp('@', line(1))
			quality = 1;
		elseif strcmp('>', line(1))
			quality = 0;
		end
	end
end

fclose(file);

