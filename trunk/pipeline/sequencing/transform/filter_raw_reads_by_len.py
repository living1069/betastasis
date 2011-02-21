#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	max_read_len = 100000;
	if len(sys.argv) == 4:
		max_read_len = int(sys.argv[3]);
	elif not 3 <= len(sys.argv) <= 4:
		print 'Bad amount of arguments. Example of proper usage:'
		print './trim_reads.py <fasta or fastq file> <min len> (<max len>)'
		sys.exit(-1)
		
	raw = open(sys.argv[1], 'r')
	min_read_len = int(sys.argv[2])
	
	while 1:
		str = fasta.readline()
		if str == '': break
		if min_read_len <= len(str)-1 <= max_read_len:
			sys.stdout.write(str)
		
if __name__ == '__main__':
	main()

