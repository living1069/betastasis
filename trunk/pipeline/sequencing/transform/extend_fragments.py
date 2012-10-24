#!/usr/bin/env python

# Author: Matti Annala <matti.annala@tut.fi>

import subprocess, sys, re

def main():
	if not len(sys.argv) == 2:
		print 'Bad amount of arguments.'
		print 'Usage: extend_fragments.py <sample.bam>'
		sys.exit(-1)
		
	bam_file = sys.argv[1]
		
	p = subprocess.Popen(
		'samtools sort -on %s %s | bedtools bamtobed -i stdin' % 
		(bam_file, bam_file), shell=True, stdout=subprocess.PIPE)
	
	prev_read = ''
	prev_chr = ''
	frag_start = -1
	frag_end = -1
	
	for line in p.stdout:
		tokens = line.split('\t')
		read_id = tokens[3][:-2]   # FIXME: Assumes /1 or /2 at the end
		if prev_read == '':
			prev_read = read_id
			prev_chr = tokens[0]
			frag_start = int(tokens[1])
			frag_end = int(tokens[2])
		elif read_id == prev_read:
			if tokens[0] == prev_chr:
				frag_start = min(frag_start, int(tokens[1]))
				frag_end = max(frag_end, int(tokens[2]))
				print '%s\t%d\t%d' % (prev_chr, frag_start, frag_end)
			prev_read = ''
			prev_chr = ''
		else:
			print '%s\t%d\t%d' % (prev_chr, frag_start, frag_end)
			prev_read = ''
			prev_chr = ''


if __name__ == '__main__':
	main()

