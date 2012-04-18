#!/usr/bin/env python

# Author: Matti Annala <matti.annala@tut.fi>

import sys


def main():
	if not len(sys.argv) == 3:
		print 'Bad amount of arguments.'
		sys.exit(-1)
		
	first = open(sys.argv[1], 'r')
	second = open(sys.argv[2], 'r')
	
	while 1:
		str_1 = first.readline()
		str_2 = second.readline()
		
		if str_1 == '' and str_2 == '': break
		
		if not len(str_1) < 1 and not str_1[0] in '#':
			if str_1[0] == '+':
				first.readline()    # Skip quality lines
			elif str_1[0] in '>@':
				prefix_1 = str_1[1:-3] if str_1[-3:-1] == '/1' else str_1[1:-1]
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
				prefix_2 = str_2[1:-3] if str_2[-3:-1] == '/1' else str_2[1:-1]
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





		
if __name__ == '__main__':
	main()

