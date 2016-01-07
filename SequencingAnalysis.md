# High throughput sequencing data analysis #

## Gene expression ##

Align transcriptome sequencing reads against the reference genome using Tophat:
```
    mkdir ../tophat_alignments
    echo *_1.fq.gz | parallel -n3 -c8 -m20 \
        'tophat2 -p8 -o ../tophat_alignments/${x%_1.fq.gz} 
         --transcriptome-index ~/tools/tophat-indexes/homo_sapiens/iGenome_37.2_NCBI
         ~/tools/bowtie2-indexes/homo_sapiens/hg19 $x ${x%_1.fq.gz}_2.fq.gz'

    cd ../tophat_alignments
    echo * | parallel -n4 \
        'samtools cat $x/{accepted_hits,unmapped}.bam > $x.bam && rm $x/{accepted_hits,unmapped}.bam'
```

Count the number of reads that aligned to the exons of each gene:
```
    sam count -b ~/organisms/homo_sapiens/ensembl_68/exons.composite.with_chr.bed *.bam
```


## Coverage over coding regions ##
```
    mkdir ../cds_coverages
    echo *.bam | parallel -n4 'sam coverage cds $x     
        ~/organisms/homo_sapiens/ensembl_68/Homo_sapiens.GRCh37.68.with_chr.gtf 
        ~/organisms/homo_sapiens/hg19.chrom.sizes > ../cds_coverages/${x%.bam}.coverage'
```