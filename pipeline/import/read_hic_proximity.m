function hic = read_hic_proximity()

hic = struct;
hic.ObsExp = cell(23, 23);
hic.RowBins = cell(23, 23);
hic.ColBins = cell(23, 23);

files = dir('.');
for k = 1:length(files)
	if files(k).name(1) == '.', continue, end
	
	tokens = regexp(files(k).name, 'HIC.*chr(.+?)_chr(.+?)_.*obsexp', 'tokens');
	if length(tokens) ~= 1, continue, end
	
	token = tokens{1};
	chrA = chromosome_sym2num(token{1});
	chrB = chromosome_sym2num(token{2});
	
	data = importdata(files(k).name);
	row_labels = data.textdata(3:end);
	col_labels = data.textdata{2};
	col_labels = textscan(col_labels, '%s', -1, 'Delimiter', '\t');
	col_labels = col_labels{1}(2:end);
	
	cols = zeros(length(col_labels), 2);
	rows = zeros(length(row_labels), 2);
	
	for k = 1:length(col_labels)
		tokens = regexp(col_labels{k}, 'chr.*:(\d+)-(\d+)', 'tokens');
		token = tokens{1};
		cols(k, :) = [str2double(token{1}) str2double(token{2})];
	end
	
	for k = 1:length(row_labels)
		tokens = regexp(row_labels{k}, 'chr.*:(\d+)-(\d+)', 'tokens');
		token = tokens{1};
		rows(k, :) = [str2double(token{1}) str2double(token{2})];
	end
	
	hic.ObsExp{chrA, chrB} = data.data;
	hic.RowBins{chrA, chrB} = rows;
	hic.ColBins{chrA, chrB} = cols;
	
	%chrA, chrB
	%whos col_labels row_labels
	%data
end


