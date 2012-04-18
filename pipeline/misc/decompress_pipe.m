function pipe = decompress_pipe(in_file)

pipe = [ptemp '.decompress_pipe'];
unix(['mkfifo ' pipe]);
unix(sprintf('gunzip -c %s > %s && rm %s &', in_file, pipe, pipe));



