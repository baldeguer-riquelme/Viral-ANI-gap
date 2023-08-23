library(tidyverse)
library(tidyr)

setwd("C:/Users/borja/OneDrive/Escritorio/Borja/Georgia_Tech/ANI_virus/fastANI/Eukaryotic_phages/Sars-Cov2")

table <- read.table("SARS-CoV-2_Metadata_NCBI.csv", header=T, sep=",")
ids_no_Ns <- read.table("ids_SARS-CoV-2_without_Ns.txt",sep="\t")

table_no_Ns <- table %>% filter(Accession %in% ids_no_Ns$V1) #Filter table to retain only sequences with no Ns

res <- table_no_Ns %>% group_by(Organization) %>% summarise(num=n())#Here I got the organization that submitted most of the genomes



#CDC sequences good quality (no Ns)
CDC_USA <- table_no_Ns %>% filter(Organization == "Centers for Disease Control and Prevention, Respiratory Viruses Branch, Division of Viral Diseases")
CDC_USA <- CDC_USA %>% mutate(WHO_variant=case_when(str_detect(Pangolin , "^B.1.1.7") ~ "Alpha",
                                                    str_detect(Pangolin , "^B.1.351") ~ "Beta",
                                                    str_detect(Pangolin , "^P.1") ~ "Gamma",
                                                    str_detect(Pangolin , "^B.1.617.2|^AY.") ~ "Delta",
                                                    str_detect(Pangolin , "^B.1.427|^B.1.429") ~ "Epsilon",
                                                    str_detect(Pangolin , "^BA.|^B.1.1.529") ~ "Omicron")) # These are the VOC, those that show a notable change in phenotype

CDC_USA_only_WHO_variants <- CDC_USA %>% drop_na(WHO_variant)
CDC_USA_only_WHO_variants %>% group_by(WHO_variant) %>% summarise(n()) 
Selected_rows_CDC_USA_only_WHO_variants <- CDC_USA_only_WHO_variants %>% group_by(WHO_variant) %>% slice_sample(n=1250)#Random extraction of 1250 genomes for each variant (Alpha, Delta, Gamma and Omicron). For Beta and Epsilon get maximum possible
accession_CDC_USA <- Selected_rows_CDC_USA_only_WHO_variants$Accession
write.table(accession_CDC_USA, file="1250_genomes_each_VOC_variant_CDC_USA.txt",sep="\t",row.names=F,col.names=F, quote=F)
write.table(Selected_rows_CDC_USA_only_WHO_variants, file="Metadata_1250_genomes_each_VOC_variant_CDC_USA.txt",sep="\t",row.names=F,col.names=F, quote=F)



#CDC sequences all genomes (with and without Ns)
CDC_USA_2 <- table %>% filter(Organization == "Centers for Disease Control and Prevention, Respiratory Viruses Branch, Division of Viral Diseases")
CDC_USA_2 <- CDC_USA_2 %>% mutate(WHO_variant=case_when(str_detect(Pangolin , "^B.1.1.7") ~ "Alpha",
                                                    str_detect(Pangolin , "^B.1.351") ~ "Beta",
                                                    str_detect(Pangolin , "^P.1") ~ "Gamma",
                                                    str_detect(Pangolin , "^B.1.617.2|^AY.") ~ "Delta",
                                                    str_detect(Pangolin , "^B.1.427|^B.1.429") ~ "Epsilon",
                                                    str_detect(Pangolin , "^BA.|^B.1.1.529") ~ "Omicron")) # These are the VOC, those that show a notable change in phenotype

CDC_USA_only_WHO_variants_all <- CDC_USA_2 %>% drop_na(WHO_variant)
CDC_USA_only_WHO_variants_all %>% group_by(WHO_variant) %>% summarise(n()) 
Selected_rows_CDC_USA_only_WHO_variants_all <- CDC_USA_only_WHO_variants_all %>% group_by(WHO_variant) %>% slice_sample(n=1250)#Random extraction of 1250 genomes for each variant (Alpha, Delta, Gamma and Omicron). For Beta and Epsilon get maximum possible
accession_CDC_USA_all <- Selected_rows_CDC_USA_only_WHO_variants_all$Accession
write.table(accession_CDC_USA_all, file="1250_genomes_each_VOC_variant_CDC_USA_with_and_without_Ns.txt",sep="\t",row.names=F,col.names=F, quote=F)
write.table(Selected_rows_CDC_USA_only_WHO_variants_all, file="Metadata_1250_genomes_each_VOC_variant_CDC_USA_with_and_without_Ns.txt",sep="\t",row.names=F,col.names=F, quote=F)




#Robert-Koch-Institute (RKI-OSEDA, Nordufer 20, 13353 Berlin) good quality sequences (no Ns)
RKI_Germany <- table_no_Ns %>% filter(Organization == "RKI-OSEDA, Nordufer 20, 13353 Berlin")
RKI_Germany <- RKI_Germany %>% mutate(WHO_variant=case_when(str_detect(Pangolin , "B.1.1.7") ~ "Alpha",
                                                    str_detect(Pangolin , "^B.1.351") ~ "Beta",
                                                    str_detect(Pangolin , "^P.1") ~ "Gamma",
                                                    str_detect(Pangolin , "^B.1.617.2|^AY.") ~ "Delta",
                                                    str_detect(Pangolin , "^B.1.427|^B.1.429") ~ "Epsilon",
                                                    str_detect(Pangolin , "^BA.|^B.1.1.529") ~ "Omicron"))# These are the VOC, those that show a notable change in phenotype

RKI_Germany_only_WHO_variants <- RKI_Germany %>% drop_na(WHO_variant)
RKI_Germany_only_WHO_variants %>% group_by(WHO_variant) %>% summarise(n()) 
Selected_rows_RKI_Germany_only_WHO_variants <- RKI_Germany_only_WHO_variants %>% group_by(WHO_variant) %>% slice_sample(n=300)#Random extraction of 300 genomes for each variant (Alpha, Delta and Omicron). For Beta, Epsilon and Gamma get maximum possible
accession_RKI_Germany <- Selected_rows_RKI_Germany_only_WHO_variants$Accession
write.table(accession_RKI_Germany, file="300_genomes_each_VOC_variant_RKI_Germany.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(Selected_rows_RKI_Germany_only_WHO_variants, file="Metadata_300_genomes_each_VOC_variant_RKI_Germany.txt",sep="\t",row.names=F,col.names=F, quote=F)




#Robert-Koch-Institute (RKI-OSEDA, Nordufer 20, 13353 Berlin) all sequences (with and without Ns)
RKI_Germany_2 <- table %>% filter(Organization == "RKI-OSEDA, Nordufer 20, 13353 Berlin")
RKI_Germany_2 <- RKI_Germany_2 %>% mutate(WHO_variant=case_when(str_detect(Pangolin , "B.1.1.7") ~ "Alpha",
                                                            str_detect(Pangolin , "^B.1.351") ~ "Beta",
                                                            str_detect(Pangolin , "^P.1") ~ "Gamma",
                                                            str_detect(Pangolin , "^B.1.617.2|^AY.") ~ "Delta",
                                                            str_detect(Pangolin , "^B.1.427|^B.1.429") ~ "Epsilon",
                                                            str_detect(Pangolin , "^BA.|^B.1.1.529") ~ "Omicron"))# These are the VOC, those that show a notable change in phenotype

RKI_Germany_only_WHO_variants_all <- RKI_Germany_2 %>% drop_na(WHO_variant)
RKI_Germany_only_WHO_variants_all %>% group_by(WHO_variant) %>% summarise(n()) 
Selected_rows_RKI_Germany_only_WHO_variants_all <- RKI_Germany_only_WHO_variants_all %>% group_by(WHO_variant) %>% slice_sample(n=300)#Random extraction of 300 genomes for each variant (Alpha, Delta and Omicron). For Beta, Epsilon and Gamma get maximum possible
accession_RKI_Germany_all <- Selected_rows_RKI_Germany_only_WHO_variants_all$Accession
write.table(accession_RKI_Germany_all, file="300_genomes_each_VOC_variant_RKI_Germany_with_and_without_Ns.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(Selected_rows_RKI_Germany_only_WHO_variants_all, file="Metadata_300_genomes_each_VOC_variant_RKI_Germany_with_and_without_Ns.txt",sep="\t",row.names=F,col.names=F, quote=F)
