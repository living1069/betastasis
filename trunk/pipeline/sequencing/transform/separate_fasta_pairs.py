#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	if len(sys.argv) != 4:
		print 'Bad amount of arguments. Example of proper usage:'
		print './separate_fasta_pairs.py <fasta_file> <5p_pair_file> <3p_pair_file>'
		sys.exit(-1)
		
	fasta = open(sys.argv[1], 'r')
	pairs_5p = open(sys.argv[2], 'w')
	pairs_3p = open(sys.argv[3], 'w')
	
	while 1:
		str = fasta.readline()
		if str == '': break
		if str[0] == '#': continue
		
		if str[0] == '>':
			header = str[1:-1]
			side_5p = (header[-2:] == '/1')
			continue
			
		if side_5p:
			pairs_5p.write('>%s\n%s' % (header, str))
		else:
			pairs_3p.write('>%s\n%s' % (header, str))
	
	pairs_5p.close()
	pairs_3p.close()
	

if __name__ == '__main__':
	main()

