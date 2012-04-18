#!/bin/bash
for num in `seq 140`; do wget -nc ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian$num.rna.gbff.gz; done
for packed in *.gz; do gunzip $packed; done

