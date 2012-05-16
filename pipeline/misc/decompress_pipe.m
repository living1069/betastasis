function pipe = decompress_pipe(in_file)

pipe = [temporary('decompress_pipe') 'pipe.fifo'];
unix(['mkfifo ' pipe]);
unix(sprintf('gunzip -c %s > %s && rm %s &', in_file, pipe, pipe));



