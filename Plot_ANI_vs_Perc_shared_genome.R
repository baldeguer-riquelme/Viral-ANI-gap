library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(data.table)

plotting <- function (tabla, species) {
  tabla['ANI_round'] <- round(tabla$ANI,1)
  tabla['Perc_round'] <- round(tabla$Perc_genome_shared,0)
  
  m<-ggplot(tabla,aes(x=ANI,y=Perc_genome_shared)) + geom_point(shape=21,fill="#F8766D", size=1.5) + 
    xlab("ANI (%)") + ylab("% shared genome") + 
    scale_x_continuous(expand=c(0,0), limits=c(94.9,100.1)) + scale_y_continuous(expand=c(0,0), limits=c(round(min(tabla$Perc_genome_shared),0)-1.1,100.5)) +
    theme_bw() + theme(text=element_text(size=12))
  
  ani_freq<-tabla %>% group_by(ANI) %>% summarise(n=n())
  
  ggplot(ani_freq, aes(x=ANI,y=n)) + geom_col(fill="#00BFC4",color="black") + 
    scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
    xlab("ANI (%)") + ylab("Frequency") +
    theme_bw() + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  
  ani_round_freq<-tabla %>% group_by(ANI_round) %>% summarise(n=n())
  
  a<-ggplot(ani_round_freq, aes(x=ANI_round,y=n)) + geom_col(fill="#00BFC4",color="black") + 
    scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0), limits=c(94.9,100.1)) +
    xlab("ANI (%)") + ylab("Frequency") +
    theme_bw() + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  
  perc_round_freq<-tabla %>% group_by(Perc_round) %>% summarise(n=n())
  
  l<-ggplot(perc_round_freq, aes(x=Perc_round,y=n)) + geom_col(fill="#00BFC4",color="black") + 
    scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0), limits=c(round(min(tabla$Perc_genome_shared),0)-1.1,100.5)) +
    xlab("% shared genome") + ylab("Frequency") + 
    theme_bw() + theme(axis.text.y  = element_blank(), axis.title.y = element_blank())
  
  pdf(paste("Plot_",species,".pdf",sep=""), width=14,height=9)
  print(ggarrange(a,NULL, m, l + coord_flip(), ncol = 2, nrow = 2, heights = c(0.25,1), widths  = c(1,0.3)))
  dev.off()
}


setwd("/storage/home/hcoda1/7/briquelme3/scratch/ANI_virus/LongReads_SRA/To_analyze")
cat(paste("[",Sys.time(),"]"," Reading ANI table\n"))
data<-read.table("ani_table_final_sp_included.txt",dec=".",sep="\t",row.names=NULL) #Table with species information. Obtained with Rscript "Species_from_ani_table"
colnames(data)<-c("Seq1","Seq2","ANI","Frags_used","Num_total_frags","Length_Seq1","Length_Seq2", "Perc_genome_shared","Species")

min_num_genomes=20 # Minimum number of genomes per species. 20 for isolated viruses, 10 for uncultured viruses (fosmid and long-reads).
cat(paste("[",Sys.time(),"]","Identifying species with more than", min_num_genomes, "genomes","\n"))

setDT(data)
Num_strain_per_sp_min10<-data.frame(Species=character(0))
for (sp in unique(data$Species)){
  #test<-subset(data, Species == sp)
  test <- data[Species == sp, ]
  column1<-as.data.frame(test$Seq1)
  column2<-as.data.frame(test$Seq2)
  colnames(column1)<-c("column")
  colnames(column2)<-c("column")
  column<-add_row(column1, column2)
  unicos<-unique(column)
  if (length(unicos$column) >= min_num_genomes){
    Num_strain_per_sp_min10$Species<-as.character(Num_strain_per_sp_min10$Species)
    Num_strain_per_sp_min10[nrow(Num_strain_per_sp_min10) + 1,] = c(sp)
    cat(paste(length(unicos$column)," strains in ", sp,"*\n",sep=""))
  } else{
    cat(paste(length(unicos$column)," strains in ", sp,"\n",sep=""))}
}
colnames(Num_strain_per_sp_min10)<-c("Species")

#Plot individual species
cat(paste("[",Sys.time(),"]","Plotting individual species plots","\n"))
npairs_min=100
npairs_max=1000000
for (sp in Num_strain_per_sp_min10$Species){ 
  plot_data<- data %>% filter(Species == sp)
  plot_data <- as.data.frame(plot_data)
  if (length(plot_data$Seq1) > npairs_min & length(plot_data$Seq1) < npairs_max) {
    cat(paste("Plotting individual species: ",sp,"\n",sep=""))
    plotting(plot_data, sp)
  } else if(length(plot_data$Seq1) > npairs_min & length(plot_data$Seq1) > npairs_max) {
    plot_data_sub <- slice_sample(plot_data, n=npairs_max)
    cat(paste("Plotting individual species: ",sp,"\n",sep=""))
    plotting(plot_data_sub, paste(sp,"subsampled",sep="_"))
  }
}

#Normalized plot. Plots all species together after randomly select a given number of pairwise comparison for each species
x=0
while (x < 2) { 
  npairs=100
  plot_data<-data.frame()
  for (sp in Num_strain_per_sp_min10$Species){ 
    sp_data<-data %>% filter(Species == sp)
    if (length(sp_data$Seq1) >= npairs){
      random_data<-sample_n(sp_data, npairs)
      plot_data<-rbind(plot_data, random_data)
    }
  }
  cat(paste("Plotting ",length(unique(plot_data$Species))," species with more than ", npairs," pairs","\n",sep=""))
  plotting(as.data.frame(plot_data), paste("All_sp_random",npairs,"pairwiserows_",x,sep=""))
  x=x+1
}


cat(paste("[",Sys.time(),"]","All done\n"))