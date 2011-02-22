#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	if len(sys.argv) != 3:
		print 'Bad amount of arguments. Example of proper usage:'
		print './trim_raw_reads.py <raw read file> <max read length>'
		sys.exit(-1)
		
	raw = open(sys.argv[1], 'r')
	max_read_len = int(sys.argv[2])

	while 1:
		str = raw.readline()
		if str == '': break
		if len(str) - 1 > max_read_len:
			sys.stdout.write(str[0:max_read_len] + '\n')
		else:
			sys.stdout.write(str)

if __name__ == '__main__':
	main()

