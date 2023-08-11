# Viral-ANI-gap

This repository contains the scripts used in the paper "A natural ANI gap that can define clusters within viral species".
Viral genomes were retrieved from NCBI and “The Actinobacteriophage Database” (https://phagesdb.org/). Accession numbers can be found in the Supplementary Table S1.


### Procedure
1. Independent viral genome sequences were placed in a folder named fastaperfile/ using the enveomics script FastA.per_file.pl (http://enve-omics.ce.gatech.edu/enveomics/docs?t=FastA.per_file.pl)

     FastA.per_file ./fastaperfile/ viral_genomes.fasta

2. A file with the path to each viral genomes was generated using:

     ls -1 fastaperfile/*.fasta > fastaperfile/list_genomes

3. fastANI in Many-to-Many mode was used to calculate ANI with a fragment length of 1kb:

      fastANI --ql ./fastaperfile/list_genomes --rl ./fastaperfile/list_genomes -o ani_table -t 32 --fragLen 1000


