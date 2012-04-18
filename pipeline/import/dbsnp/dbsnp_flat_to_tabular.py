
import os
import re
import gzip

snp_files = [f for f in os.listdir('.') if re.match('.*ch.+.flat.gz', f)]
for filename in snp_files:
	f = gzip.open(filename)
	
	print 'CHROMOSOME\tCHROM_POS\tSNP_ID\tALLELES\tHETEROZYGOSITY'
	
	for line in f:
		m = re.search(r'^(rs\d+)', line)
		if m:
			snp_id = m.group(1)
		#if re.search('^ss.+1000GENOMES'):
		
		m = re.search(r'^SNP \| alleles=\'(.+?)\' \| het=(.+?) \|', line)
		if m:
			alleles = m.group(1)
			heterozygosity = m.group(2)
		
		m = re.search(r'^CTG \| .+GRCh37.+ \| chr=(.+?) \| chr-pos=(\d+)', line)
		if m:
			chromosome = m.group(1)
			chr_pos = m.group(2)
			print 'chr%s\t%s\t%s\t%s\t%s' % (
				chromosome, chr_pos, snp_id, alleles, heterozygosity)
	
	f.close()



