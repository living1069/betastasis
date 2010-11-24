#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def split_fasta_read(fasta, tag_length, color):
	str = ''
	header = ''
	
	while 1:
		header = fasta.readline()
		if header == '': return 1
		if header[0] == '>' or header[0] == '@':
			break
	
	while 1:
		str = fasta.readline()
		if str == '':
			return 1
		if str[0] == '#': continue
			
		if len(str) < 2 * tag_length + 5:
			break
		
		sys.stdout.write(header[:-1] + '~1\n')
		sys.stdout.write(str[0:tag_length+1] + '\n')
		sys.stdout.write(header[:-1] + '~2\n')
		if color:
			sys.stdout.write('T' + str[-tag_length-1:])
		else:
			sys.stdout.write(str[-tag_length:])
		break
	
	if header[0] == '>': return 0
	
	print 'Tag pair splitting does not work for FASTQ files at the moment.'
	return 1




def main():
	if len(sys.argv) != 4:
		print 'Bad amount of arguments. Example of proper usage:'
		print './split_reads_into_tag_pairs.py <fasta_file> <tag_len> [nucleotide|color]'
		sys.exit(-1)
		
	color = 0
	
	if sys.argv[3] == 'nucleotide':
		color = 0
	elif sys.argv[3] == 'color':
		color = 1
	else:
		print 'Invalid read type argument.'
		sys.exit(-1)
	
	fasta = open(sys.argv[1], 'r')
	tag_length = int(sys.argv[2])

	while 1:
		if split_fasta_read(fasta, tag_length, color):
			break
		
if __name__ == '__main__':
	main()

