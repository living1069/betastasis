#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	if len(sys.argv) != 5:
		print 'Bad amount of arguments. Example of proper usage:'
		print './split_reads_into_tag_pairs.py <raw_file> <tag_len> <5p_tag_file> <3p_tag_file>'
		sys.exit(-1)
		
	color = -1
	
	raw = open(sys.argv[1], 'r')
	tag_length = int(sys.argv[2])
	5p_tags = open(sys.argv[3], 'w')
	3p_tags = open(sys.argv[4], 'w')

	while 1:
		str = raw.readline()
		if str == '': break
		if str[0] == '#': continue
		if len(str) < 2 * tag_length + 5: continue
		
		if color == -1:
			color = (str[1] in '0123.')
		
		if color:
			5p_tags.write(str[0:tag_length+1] + '\n')
			3p_tags.write('T' + str[-tag_length-1:])
		else:
			5p_tags.write(str[0:tag_length] + '\n')
			3p_tags.write(str[-tag_length:])
	
	5p_tags.close()
	3p_tags.close()


if __name__ == '__main__':
	main()

