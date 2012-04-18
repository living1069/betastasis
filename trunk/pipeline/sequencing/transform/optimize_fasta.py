#!/usr/bin/env python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def fasta_compact_single(fasta_file):
	read_id = 0
	fasta = open(fasta_file, 'r')
	while 1:
		str = fasta.readline()
		if str == '': break
		if len(str) < 1 or str[0] in '#>@': continue
		if str[0] == '+':
			fasta.readline()
			continue
		read_id += 1
		sys.stdout.write('>%d\n' % read_id)
		sys.stdout.write(str)


def fasta_compact_paired(fastq_1_file, fastq_2_file,
	fasta_1_file, fasta_2_file):
	
	first = open(fastq_1_file, 'r')
	second = open(fastq_2_file, 'r')
	
	first_raw = open(fasta_1_file, 'w')
	second_raw = open(fasta_2_file, 'w')
	
	first_prefixes = {}
	second_prefixes = {}
	prefix_1 = ''
	prefix_2 = ''
	
	read_id = 0
	
	while 1:
		str_1 = first.readline()
		str_2 = second.readline()
		
		if str_1 == '' and str_2 == '': break
		
		if not len(str_1) < 1 and not str_1[0] in '#':
			if str_1[0] == '+':
				first.readline()    # Skip quality lines
			elif str_1[0] in '>@':
				if str_1[-3:-1] == '/1': prefix_1 = str_1[1:-3]
				else: prefix_1 = str_1[1:-1]
			else:
				pair = second_prefixes.pop(prefix_1, -1)
				if pair != -1:
					read_id += 1
					first_raw.write('>%d\n' % read_id)
					first_raw.write(str_1)
					second_raw.write('>%d\n' % read_id)
					second_raw.write(pair)
				else:
					first_prefixes[prefix_1] = str_1
		
		if not len(str_2) < 1 and not str_2[0] in '#':
			if str_2[0] == '+':
				second.readline()     # Skip quality lines
			elif str_2[0] in '>@':
				if str_2[-3:-1] == '/2': prefix_2 = str_2[1:-3]
				else: prefix_2 = str_2[1:-1]
			else:
				pair = first_prefixes.pop(prefix_2, -1)
				if pair != -1:
					read_id += 1
					first_raw.write('>%d/1\n' % read_id)
					first_raw.write(pair)
					second_raw.write('>%d/2\n' % read_id)
					second_raw.write(str_2)
				else:
					second_prefixes[prefix_2] = str_2
	
	if len(first_prefixes) > 0:
		print 'Found %d first reads without a pair.' % len(first_prefixes)
	if len(second_prefixes) > 0:
		print 'Found %d second reads without a pair.' % len(second_prefixes)
	
	first_raw.close()
	second_raw.close()


def main():
	if not len(sys.argv) in [2, 5]:
		print 'Bad amount of arguments.'
		sys.exit(-1)
		
	if len(sys.argv) == 2:
		fasta_compact_single(sys.argv[1])
	else:
		fasta_compact_paired(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

		
if __name__ == '__main__':
	main()

