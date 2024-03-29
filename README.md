# Viral-ANI-gap

This repository contains the scripts and an extended version of the procedure used in the paper "A natural ANI gap that can define clusters within viral species".
Viral genomes were retrieved from NCBI and “The Actinobacteriophage Database” (https://phagesdb.org/). Accession numbers can be found in the Supplementary Table S1.


### General procedure
1. Independent viral genome sequences were placed in a folder named fastaperfile/ using the enveomics script FastA.per_file.pl (http://enve-omics.ce.gatech.edu/enveomics/docs?t=FastA.per_file.pl)

     ```FastA.per_file ./fastaperfile/ viral_genomes.fasta```

2. A file with the path to each viral genomes was generated using:

     ```ls -1 fastaperfile/*.fasta > fastaperfile/list_genomes```

3. fastANI in Many-to-Many mode was used to calculate ANI with a fragment length of 1kb:

     ```fastANI --ql ./fastaperfile/list_genomes --rl ./fastaperfile/list_genomes -o ani_table -t 12 --fragLen 1000```

4. Genome length was calculated using the enveomics script FastA.length.pl (http://enve-omics.ce.gatech.edu/enveomics/docs?t=FastA.length.pl):

```FastA.length.pl viral_genomes.fasta > Length_viruses.txt```

5. Output generated on step 3 by fastANI (ani_table) and on step 4 by FastA.length (Length_viruses.txt) were then used as input for the script "Species_from_ani_table.R" which does the following steps:
   
   5.1 Select pairs above 95% ANI (threshold for viral species).
   
   5.2 Cluster the sequences sharing more than this value to identify species.
   
   5.3 Remove self-hits.

   5.4 Keep only species with more than 20 or 10 genomes for cultured and uncultured viral species, respectively.
   
   5.5 Write output named "ani_table_final_sp_included.txt"

```Rscript Species_from_ani_table.R```

6. The script "Plot_ANI_vs_Perc_shared_genome.R" was then used to represent the ANI vs percentage of shared genome. The script does the following steps:
   
   6.1 Plot individual species
   
   6.2 Plot all species together after normalizing the data to 150 pairs per species

```Rscript Plot_ANI_vs_Perc_shared_genome.R```
   
![alt text](https://github.com/baldeguer-riquelme/Viral-ANI-gap/blob/main/.figure/Figure1.png)



Plots for each species can be found in folders "Plots_prokaryotic_viruses", "Plots_eukaryotic_viruses", "Plots_fosmid_viruses" and "Plots_long_read_viruses"


### Bootstrapped subsampling and peak detection
The script "Bootstrap_peak_detector.py" uses the "ani_table_final_sp_included.txt" file (outputted by the "Species_from_ani_table.R" script) as input and does the following steps:

1.1 It performs a random subsampling of the data to get the same number of pairs per species and identify peaks and valleys. This process is repeated 1,000 times, the results are pooled and then plotted.
     
1.2 The output is a pdf called "Bootstrap_plot.pdf", with two plots: at the bottom, the ANI values distribution and the positions of the peaks and valleys; at the top, an histogram displaying the number of times each peak and valley was detected.

```python Bootstrap_peak_detector.py```


### Species classification into 4 groups
The script "Peak_detector.py" uses the "ani_table_final_sp_included.txt" file (outputted by the "Species_from_ani_table.R" script) as input and does the following steps:

1.1 For each species, the script extracts and smooth the data using the gaussian_kde function of the scipy package.
     
1.2 Then, it assess normality and unimodality using the D’Agostino and Pearson’s test and Hartigan’s dip test, respectively.
     
1.3 The smooth data is used to identify peaks and valleys with the function find_peaks() of the scipy package.
     
1.4 These results are integrated and the species is classified into one of the four groups described in the paper. The result is printed to stdout.
     
1.5 Finally, it produces a pdf named "Peak_Valley_per_species.pdf" with plots for all species which includes the original data histogram, smoothed distribution, peaks and valleys detected, p-values of normality and unimodality tests (top-right) and the group the species belong to (top-left, below species name).

```python Peak_detector.py```

### SARS-CoV-2
SARS-CoV-2 genomes were downloaded from NCBI as well as their associated metadata. Since the amount of genomes is too big, a subset of sequences were selected using the script "SARS-CoV-2_select_genomes.R". This script was used to:

1. identify the main submitting organizations, that were the Center for Disease Control (CDC) in the USA and the Robert-Koch-Institute in Germany.
     
2. randomly select a given number of genomes for each Variant of Concern (Alpha, Beta, Delta, Gamma, Epsilon and Omicron) and submitter.
     
3. write a file with the randomly selected genome ids.
     
4. write a file with the metadata associated with the randomly selected genomes.

This procedure of random selection was used to subset all genomes (with or without Ns) as well as the high-quality genomes (no Ns, identified previously by FastA.filterN.pl from the enveomics collection: http://enve-omics.ce.gatech.edu/enveomics/).


Selected sequences were then extracted from the complete data set using the "FastA.filter.pl" script from enveomics and then, the general procedure explained above was followed. However, in this case, on step 5, instead of the script "Plot_ANI_vs_Perc_shared_genome.R", the script "Processing_script_ncbi_Sars-Cov2_CDC_USA.R" was employed to plot histograms with and without the variant data and to calculate the percentage of comparisons between the same or different genomovar above and below the ANI gap.

![alt text](https://github.com/baldeguer-riquelme/Viral-ANI-gap/blob/main/.figure/Histogram_Sars-CoV-2_variants.png)
