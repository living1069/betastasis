#!/usr/bin/python

import os
import re

def parse_chr(source, target):
	file = open(source, 'r')
	writing = False
	for line in file:
		if line[0] == '>':
			if writing: writing = False
				
			m = re.search('\|NC_(\d+)\.\d+\|', line)
			if not m: continue
				
			id = int(m.group(1))
			if id >= 1 and id <= 24:
				sym = ('chr%d' % id)
				if id == 23: sym = 'chrX'
				if id == 24: sym = 'chrY'
					
				print 'Found chromosome %s.' % sym
				
				writing = True
				target.write('>%s\n' % sym)
		elif writing:
			target.write(line)
	
	file.close()

fna_in_files = [ f for f in os.listdir('.')
	if re.match('vertebrate.*\.genomic\.fna', f) ]

target = open('genome.fa', 'w')

for filename in fna_in_files:
	parse_chr(filename, target)

