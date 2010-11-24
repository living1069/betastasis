#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def trim_qual(fastq):
	str = ''
	header = ''
	
	while 1:
		header = fastq.readline()
		if header == '':
			return 1
		if header[0] != '#':
			out = '>' + header[1:]
			sys.stdout.write(out)
			break
	
	while 1:
		str = fastq.readline()
		if str == '':
			return 1
		if str[0] != '#':
			sys.stdout.write(str)
			break
	
	while 1:
		str = fastq.readline()
		if str == '':
			return 1
		if str[0] == '#':
			continue
		break
	
	while 1:
		str = fastq.readline()
		if str == '':
			return 1
		if str[0] == '#':
			continue
		break
	
	return 0


def main():
	if len(sys.argv) != 2:
		print 'Bad amount of arguments.'
		sys.exit(-1)
		
	fastq = open(sys.argv[1], 'r')

	while 1:
		if trim_qual(fastq):
			break
		
if __name__ == '__main__':
	main()

