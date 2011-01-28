#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	if len(sys.argv) != 2:
		print 'Bad amount of arguments. Example of proper usage:'
		print './raw_to_fasta.py <raw read file>'
		sys.exit(-1)
		
	raw = open(sys.argv[1], 'r')

	str = ''
	read_num = 1
	
	while 1:
		str = raw.readline()
		if str == '': break
		if len(str) < 2: continue
			
		sys.stdout.write('>%d\n' % read_num)
		sys.stdout.write(str)
		read_num += 1

		
if __name__ == '__main__':
	main()

