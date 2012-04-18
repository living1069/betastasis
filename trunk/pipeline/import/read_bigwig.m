function [chr, start, stop, val] = read_bigwig(bigwig_file)

assume_ordered = true;

pipe = tempname();
%sprintf('mkfifo %s && bigWigToWig %s >(wig2bed - > %s) &', pipe, bigwig_file, pipe)
unix(sprintf('mkfifo %s && bigWigToWig %s >(wig2bed - > %s) &', ...
	pipe, bigwig_file, pipe));

pause(1);  % Wait for mkfifo to finish, otherwise function may fail.

fid = fopen(pipe);

if assume_ordered
	data = textscan(fid, '%*s %d %d %*s %f', 'Delimiter', '\t', ...
		'ReturnOnError', 0);
	chr = nan(length(data{1}), 1);
	start = data{1};
	stop = data{2};
	val = data{3};
else
	data = textscan(fid, '%s %d %d %*s %f', 'Delimiter', '\t', ...
		'ReturnOnError', 0);
	chr = chromosome_sym2num(data{1});
	start = data{2};
	stop = data{3};
	val = data{4};
end
fclose(fid);

delete(pipe);

