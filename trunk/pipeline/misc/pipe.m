function pipe = pipe(command, pipe)

unix(sprintf('rm %s 2> /dev/null; mkfifo %s', pipe, pipe));
unix(sprintf('%s >> %s &', command, pipe));
pause(1);

