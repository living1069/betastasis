#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def merge_read(fasta, max_read_len):
	str = ''
	header = ''
	
	while 1:
		header = fasta.readline()
		if header == '': return 1
		if header[0] == '>' or header[0] == '@':
			sys.stdout.write(header)
			break
	
	while 1:
		str = fasta.readline()
		if str == '': return 1
		if str[0] != '#':
			if len(str) - 1 > max_read_len:
				sys.stdout.write(str[0:max_read_len] + '\n')
			else:
				sys.stdout.write(str)
			break
		
	return 0
	
	# FIXME: This old code only works on SOLiD FASTQ files.
	if header[0] == '>': return 0
	
	while 1:
		str = fasta.readline()
		if str == '': return 1
		if str[0] == '+':
			sys.stdout.write(str)
			break
	
	while 1:
		str = fasta.readline()
		if str == '': return 1
		if str[0] == '#': continue
			
		pos = 0
		for k in range(1, 18):
			pos = str.find(' ', pos) + 1
		pos = str.find(' ', pos)
		sys.stdout.write(str[0:pos] + '\n')
		break
	
	return 0


def main():
	if len(sys.argv) != 3:
		print 'Bad amount of arguments. Example of proper usage:'
		print './trim_reads.py <fasta or fastq file> <max read length>'
		sys.exit(-1)
		
	fasta = open(sys.argv[1], 'r')
	max_read_len = int(sys.argv[2])

	while 1:
		if merge_read(fasta, max_read_len):
			break
		
if __name__ == '__main__':
	main()

