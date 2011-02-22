#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def filter_read(fasta, min_read_len, max_read_len):
	while 1:
		header = fasta.readline()
		if header == '': return 1
		if header[0] == '>':
			break
	
	while 1:
		str = fasta.readline()
		if str == '': return 1
		if str[0] == '#': continue
		read_len = len(str) - 1
		if read_len >= min_read_len and read_len <= max_read_len:
			sys.stdout.write(header)
			sys.stdout.write(str)
		break
	
	return 0
	


def main():
	max_read_len = 100000;
	if len(sys.argv) == 4:
		max_read_len = int(sys.argv[3]);
	elif not 3 <= len(sys.argv) <= 4:
		print 'Bad amount of arguments. Example of proper usage:'
		print './trim_reads.py <fasta or fastq file> <min len> (<max len>)'
		sys.exit(-1)
		
	fasta = open(sys.argv[1], 'r')
	min_read_len = int(sys.argv[2])
	
	while 1:
		if filter_read(fasta, min_read_len, max_read_len):
			break
		
if __name__ == '__main__':
	main()

