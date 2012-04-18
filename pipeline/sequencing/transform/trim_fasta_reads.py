#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	if len(sys.argv) not in [3, 4]:
		print 'Bad amount of arguments. Example of proper usage:'
		print './trim_fasta_reads.py <fasta or fastq file> <max read length> [suffix]'
		sys.exit(-1)
		
	fasta = open(sys.argv[1], 'r')
	trim_len = int(sys.argv[2])
	suffix = sys.argv[3] if len(sys.argv) == 4 else ''
	
	header = ''
	color = -1
	
	trim_3p = False
	if trim_len < 0:
		trim_3p = True
		trim_len = -trim_len
	
	while 1:
		str = fasta.readline()
		if str == '': break
		if str[0] == '#': continue
		if str[0] == '>':
			header = str[:-1]
			continue
		
		if color == -1:
			color = (str[1] in '0123.')

		sys.stdout.write(header + suffix + '\n')
		
		if len(str) - 1 <= trim_len:
			sys.stdout.write(str)
		elif trim_3p:
			seq_prefix = 'T' if color else ''
			sys.stdout.write(seq_prefix + str[-trim_len:-1] + '\n')
		else:
			sys.stdout.write(str[0:trim_len] + '\n')

		
if __name__ == '__main__':
	main()

