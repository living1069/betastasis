#!/usr/bin/python

# Author: Matti Annala <matti.annala@tut.fi>
#
# This script should be executed inside a directory containing a full RefSeq
# release with all vertebrate_mammalian_* files. The script will collect all
# transcript and genomic annotations for the species Homo sapiens.
# 
# The script will write the collected annotations into the following files:
#  * refseq_human.rna.gbff
#  * refseq_human.genomic.gbff
#

import os
import re

def merge_gbff(target, source):
	print 'Collecting annotations from file \'%s\'...' % source
	file = open(source, 'r')
	transcript = ''
	selected = False
	for line in file:
		transcript += line

		if line == '//\n':
			if selected:
				target.write(transcript)
			transcript = ''
			selected = False
			continue
		
		if re.match('.*ORGANISM.*Homo sapiens', line):
			selected = True
			continue
			
	file.close()

rna_gbff_files = [ f for f in os.listdir('.')
	if re.match('vertebrate.*\.rna\.gbff', f) ]
genomic_gbff_files = [ f for f in os.listdir('.')
	if re.match('vertebrate.*\.genomic\.gbff', f) ]

rna_gbff_human = open('refseq_human.rna.gbff', 'w')
print 'Compiling transcript annotations in GBFF format:'
for filename in rna_gbff_files:
	merge_gbff(rna_gbff_human, filename)
rna_gbff_human.close()

genomic_gbff_human = open('refseq_human.genomic.gbff', 'w')
print 'Compiling genomic annotations in GBFF format:'
for filename in genomic_gbff_files:
	merge_gbff(genomic_gbff_human, filename)
genomic_gbff_human.close()
