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
	
	read_id = 0

	while 1:
		str = fasta.readline()
		if str == '': break
		if str[0] in '#>': continue
		
		if len(str) < 2 * tag_length + 5: continue
		
		if color == -1:
			color = (str[1] in '0123.')
		
		read_id += 1
		
		if color:
			tags_5p.write('>%d/1\n' % read_id)
			tags_5p.write(str[0:tag_length+1] + '\n')
			tags_3p.write('>%d/2\n' % read_id)
			tags_3p.write('T' + str[-tag_length-1:])
		else:
			tags_5p.write('>%d/1\n' % read_id)
			tags_5p.write(str[0:tag_length] + '\n')
			tags_3p.write('>%d/2\n' % read_id)
			tags_3p.write(str[-tag_length-1:])
	
	tags_5p.close()
	tags_3p.close()
	

if __name__ == '__main__':
	main()

