#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def main():
	if not len(sys.argv) in [2, 5]:
		print 'Bad amount of arguments.'
		sys.exit(-1)
		
	if len(sys.argv) == 2:
		fastq = open(sys.argv[1], 'r')
		while 1:
			str = fastq.readline()
			if str == '': break
			if len(str) < 1 or str[0] in '#>@': continue
			if str[0] == '+':
				fast.readline()
				continue
			sys.stdout.write(str)
	else:
		first = open(sys.argv[1], 'r')
		second = open(sys.argv[2], 'r')
		
		first_raw = open(sys.argv[3], 'w')
		second_raw = open(sys.argv[4], 'w')
		
		first_prefixes = {}
		second_prefixes = {}
		prefix_1 = ''
		prefix_2 = ''
		
		while 1:
			str_1 = first.readline()
			str_2 = second.readline()
			
			if str_1 == '' and str_2 == '': break
			
			if not len(str_1) < 1 and not str_1[0] in '#':
				if str_1[0] == '+':
					first.readline()    # Skip quality lines
				elif str_1[0] in '>@':
					if not str_1[-3:-1] == '/1':
						print 'ERROR: Missing pair tag /1 (%s)' % str_1
						sys.exit(-1)
					prefix_1 = str_1[1:-3]
				else:
					pair = second_prefixes.pop(prefix_1, -1)
					if pair != -1:
						first_raw.write(str_1)
						second_raw.write(pair)
					else:
						first_prefixes[prefix_1] = str_1
			
			if not len(str_2) < 1 and not str_2[0] in '#':
				if str_2[0] == '+':
					second.readline()     # Skip quality lines
				elif str_2[0] in '>@':
					if not str_2[-3:-1] == '/2':
						print 'ERROR: Missing pair tag /2 (%s).' % str_2
						sys.exit(-1)
					prefix_2 = str_2[1:-3]
				else:
					pair = first_prefixes.pop(prefix_2, -1)
					if pair != -1:
						first_raw.write(pair)
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

