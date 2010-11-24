function len_dist = seq_read_lengths(reads, image)

[color, quality] = seq_read_type(reads);

file = fopen(reads, 'r');

len_dist = zeros(1, 500);
expect_header = 1;

while 1
	line = fgetl(file);
	if line == -1, break, end
		
	if length(line) == 0 || strcmp('#', line(1)), continue, end

	if expect_header
		if quality && line(1) == '@'
			expect_header = 0;
		elseif ~quality && line(1) == '>'
			expect_header = 0;
		end
	else
		if regexpi(line, '^[tcgarykmswbdhvn]+$')
			len = length(line);
			len_dist(len) = len_dist(len) + 1;
			expect_header = 1;
		elseif regexpi(line, '^[ACTG][0123\.]+$')
			len = length(line) - 1;
			len_dist(len) = len_dist(len) + 1;
			expect_header = 1;
		end
	end
end

fclose(file);

if nargin == 2
	fprintf(1, 'Rendering a histogram of the read offset distribution...\n');
	figure;
	bar(len_dist); xlim([0 (max(find(len_dist ~= 0)) + 1)]);
	saveas(gcf, image);
end

