#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def merge_read(fasta, qual):
	str = ''
	header = ''
	
	while 1:
		header = fasta.readline()
		if header == '':
			return 1
		if header[0] != '#':
			out = '@' + header[1:]
			sys.stdout.write(out)
			break
	
	while 1:
		str = fasta.readline()
		if str == '':
			return 1
		if str[0] != '#':
			sys.stdout.write(str)
			break
	
	while 1:
		str = qual.readline()
		if str == '':
			return 1
		if str[0] == '#':
			continue
		if str != header:
			sys.stderr.write('Header mismatch found!\n')
			return 1
		out = '+' + str[1:]
		sys.stdout.write(out)
		break
	
	while 1:
		str = qual.readline()
		if str == '':
			return 1
		if str[0] == '#':
			continue
		
		sys.stdout.write(str)
		break
	
	return 0


def main():
	if len(sys.argv) != 3:
		print 'Bad amount of arguments.'
		sys.exit(-1)
		
	fasta = open(sys.argv[1], 'r')
	qual = open(sys.argv[2], 'r')

	while 1:
		if merge_read(fasta, qual):
			break
		
if __name__ == '__main__':
	main()

