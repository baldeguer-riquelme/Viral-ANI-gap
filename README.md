# Viral-ANI-gap

This repository contains the scripts and an extended version of the procedure used in the paper "A natural ANI gap that can define clusters within viral species".
Viral genomes were retrieved from NCBI and “The Actinobacteriophage Database” (https://phagesdb.org/). Accession numbers can be found in the Supplementary Table S1.


### General procedure
1. Independent viral genome sequences were placed in a folder named fastaperfile/ using the enveomics script FastA.per_file.pl (http://enve-omics.ce.gatech.edu/enveomics/docs?t=FastA.per_file.pl)

     FastA.per_file ./fastaperfile/ viral_genomes.fasta


2. A file with the path to each viral genomes was generated using:

     ls -1 fastaperfile/*.fasta > fastaperfile/list_genomes


3. fastANI in Many-to-Many mode was used to calculate ANI with a fragment length of 1kb:

      fastANI --ql ./fastaperfile/list_genomes --rl ./fastaperfile/list_genomes -o ani_table -t 32 --fragLen 1000


4. Output generated by fastANI (ani_table) was then used as input for the script "Species_from_ani_table.R" which does the following steps:
   
   4.1 Select pairs above 95% ANI (threshold for viral species).
   
   4.2 Cluster the sequences sharing more than this value to identify species.
   
   4.3 Remove self-hits.
   
   4.4 Write output named ani_table_final_sp_included.txt
   


5. The script "Plot_ANI_vs_Perc_shared_genome.R") was then used to represent the ANI vs percentage of shared genome for species with more than 20 or 10 genomes for cultured and uncultured viral species, respectively. The script does the following steps:
   
   5.1 Filter species with insufficient number of genomes
   
   5.2 Plot individual species
   
   5.3 Plot all species together after normalizing the data to 150 pairs per species
   
![alt text](https://github.com/baldeguer-riquelme/Viral-ANI-gap/blob/main/.figure/Figure1.tiff)

Plots for each species can be found in folders "Plots_prokaryotic_viruses", "Plots_eukaryotic_viruses", "Plots_fosmid_viruses" and "Plots_long_read_viruses"
