function pipe = compress_pipe(out_file)

pipe = [ptemp '.compress_pipe'];
unix(['mkfifo ' pipe]);
unix(sprintf('cat %s | gzip -c > %s && rm %s &', pipe, out_file, pipe));


