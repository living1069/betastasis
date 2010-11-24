#!/bin/bash
for num in `seq 56`; do wget -nc ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian$num.rna.gbff.gz; done
for num in `seq 56`; do wget -nc ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian$num.genomic.gbff.gz; done
wget -nc ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian18.genomic.fna.gz;
wget -nc ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/H_sapiens/Assembled_chromosomes/hs_ref_GRCh37_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.agp.gz
for packed in *.gz; do gunzip $packed; done
