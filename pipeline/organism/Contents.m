% ORGANISM
%
%    The pipeline provides biological information about a number of different
%    organisms. This knowledge is stored in versioned data structures known as 
%    "organism builds". The pipeline always has one organism build which
%    is accessible through the global variable "organism". When the pipeline is
%    started, the default organism build is loaded, but the user can select
%    another organism or version using the function SELECT_ORGANISM. Available
%    organism builds can be listed using LIST_ORGANISMS.
%    
%    The currently selected (active) organism build is used by the pipeline
%    when running analyses. For example, gene expression probesets you construct
%    are created for the currently selected organism build. All datasets in the
%    pipeline are also associated with a particular organism build. This is done
%    to prevent situations where datasets based on different biological
%    assumptions were to be used together, leading to anomalous results.
%
%    The information provided by organism builds can be categorized as follows:
%    - Genetic information (DNA/RNA annotations, relationships, variations)
%    - Biological networks (gene regulatory networks, pathways, metabolism)
%    - Phenotype <-> genotype associations
%    - A variety of different ontologies
%
%    When an organism build is first loaded by the pipeline, not all of the
%    information is immediately loaded into memory. Often an organism build
%    can contain gigabytes of data, most of which is only rarely used in
%    analyses. The pipeline instead loads the organism data into memory only
%    when first accessed by the user or some analysis function. The data then
%    stays cached in memory indefinitely. This means that when you for instance
%    first access the SNP information for chromosome 1, it will take a couple
%    of seconds before the data is loaded from the disk onto memory. Subsequent
%    accesses to these SNPs will then be served quickly, using cached data.
%
%    GENOME/TRANSCRIPTOME ANNOTATIONS:
%    ---------------------------------
%    An organism's genes are described by the structure organism.Genes. For each
%    gene the structure provides the gene name and Entrez ID. The number of
%    different RNAs (splice variants) transcribed from a gene can be seen in
%    field organism.Genes.TranscriptCount. The field organism.Genes.Transcripts
%    contains indices into the transcript annotations in organism.Transcripts.
%
%    The fields in organism.Transcripts provide a name and sequence for each
%    RNA transcript. They also describe the exons of which the transcript
%    is composed of, and the position of the coding region (CDS) relative
%    to the transcripts beginning. The field organism.Transcripts.Gene links
%    each transcript to the gene from which it is transcribed, by providing
%    an index to the structure organism.Genes.
%    
