#!/usr/bin/env python

# Author: Matti Annala <matti.annala@tut.fi>

import sys
import os
import re

def main():
	if len(sys.argv) < 2:
		print 'You need to specify one or more .bam files. Example:'
		print './bam_to_fasta.py <bam file>'
		sys.exit(-1)
	
	tmp_file = '/worktmp/tmp/askfldjsal.sam'
	
	unmapped_regexp = re.compile(r'^(\S+)\t4\t.+CS:Z:([ACGT][0-9.]+)\t')
	mapped_regexp = re.compile(r'^(\S+)\t.+CS:Z:([ACGT][0-9.]+)\t')
	
	read_ids = set()
	
	for bam_file in sys.argv[1:]:
		os.system('/worktmp/pipeline/tools/samtools/samtools view ' +
			bam_file + ' > ' + tmp_file)
		
		sam = open(tmp_file, 'r')
		
		while 1:
			line = sam.readline()
			if line == '': break
			
			tokens = line.split('\t')
			if tokens[1] != '4':
				if tokens[0] in read_ids: continue
				read_ids.add(tokens[0]) 
				
			sys.stdout.write('%s\n' % tokens[9])
				
			#m = unmapped_regexp.search(line)
			#if m:
			#	sys.stdout.write('>%s\n' % m.group(1))
			#	sys.stdout.write('%s\n' % m.group(2))
			#	continue
			
			#m = mapped_regexp.search(line)
			#if m:
			#	if m.group(1) in read_ids: continue
			#	read_ids.add(m.group(1)) 
			#	sys.stdout.write('>%s\n' % m.group(1))
			#	sys.stdout.write('%s\n' % m.group(2))
			#	continue
		
		sam.close()
		

		
if __name__ == '__main__':
	main()

