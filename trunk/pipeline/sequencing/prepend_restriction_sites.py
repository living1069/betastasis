#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>

import sys

def prepend_to_fastq(fastq, restriction_site):
	str = ''
	header = ''
	
	while 1:
		header = fastq.readline()
		if header == '':
			return 1
		if header[0] == '@':
			sys.stdout.write('>' + header[1:])
			break
	
	while 1:
		str = fastq.readline()
		if str == '':
			return 1
		if str[0] == '#': continue
		
		sys.stdout.write(restriction_site + str)
		break
	
	while 1:
		str = fastq.readline()
		if str == '':
			return 1
		if str[0] == '+':
			break
	
	while 1:
		str = fastq.readline()
		if str == '':
			return 1
		if str[0] == '#': continue
		break
	
	return 0


def main():
	if len(sys.argv) != 3:
		print 'Bad amount of arguments. Example of proper usage:'
		print './prepend_restriction_sites.py <fastq_file> <restriction_site_seq>'
		sys.exit(-1)
	
	fastq = open(sys.argv[1], 'r')
	restriction_site = sys.argv[2]

	while 1:
		if prepend_to_fastq(fastq, restriction_site):
			break
		
if __name__ == '__main__':
	main()

