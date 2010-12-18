#!/usr/bin/python

import sys

def main():
	if len(sys.argv) != 2:
		print 'Bad amount of arguments. Example of proper usage:'
		print './split_fasta.py <fasta file>'
		sys.exit(-1)
		
	fasta = open(sys.argv[1], 'r')
	
	seq_file = ''
	
	while 1:
		line = fasta.readline()
		if line == '': return 0
		if line[0] == '>':
			if seq_file:
				seq_file.close()
				seq_file = ''
			seq_file = open(line[1:-1] + '.seq', 'w')
		else:
			line.replace(' ', '')
			seq_file.write(line[0:-1])
	
	seq_file.close()
		
if __name__ == '__main__':
	main()

