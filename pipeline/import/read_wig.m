function [chr, pos, val] = read_wig(wig_file)

pipe = tempname();
unix(sprintf('mkfifo %s && swiss wig2bed %s > %s &', pipe, wig_file, pipe));

pause(1);  % Wait for mkfifo to finish, otherwise function may fail.

fid = fopen(pipe);

data = textscan(fid, '%s %d %*d %*s %f', 'Delimiter', '\t', 'ReturnOnError', 0);
chr = chromosome_sym2num(data{1});
pos = data{2};
val = data{3};
fclose(fid);

delete(pipe);

