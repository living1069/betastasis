#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	if len(sys.argv) != 5:
		print 'Bad amount of arguments. Example of proper usage:'
		print './split_fasta_reads.py <fasta_file> <tag_len> <5p_tag_file> <3p_tag_file>'
		sys.exit(-1)
		
	color = -1
	
	fasta = open(sys.argv[1], 'r')
	tag_length = int(sys.argv[2])
	tags_5p = open(sys.argv[3], 'w')
	tags_3p = open(sys.argv[4], 'w')
	
	while 1:
		str = fasta.readline()
		if str == '': break
		if str[0] == '#': continue
		if str[0] == '>':
			header = str[:-1]
			continue
		
		if len(str) < 2 * tag_length + 2: continue
		
		if color == -1:
			color = (str[1] in '0123.')
		
		if color:
			tags_5p.write('%s/1\n%s\n' % (header, str[0:tag_length+1]))
			tags_3p.write('%s/2\nT%s\n' % (header, str[-tag_length-1:]))
		else:
			tags_5p.write('%s/1\n%s\n' % (header, str[0:tag_length]))
			tags_3p.write('%s/2\n%s\n' % (header, str[-tag_length-1:]))
	
	tags_5p.close()
	tags_3p.close()
	

if __name__ == '__main__':
	main()

