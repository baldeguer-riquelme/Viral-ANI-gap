library(network)
library(sna)
library(reshape2)
library(dplyr)
library(tidyr)
library(igraph)
library(ggpubr)
library(data.table)

# setwd() if needed

cat(paste("[",Sys.time(),"]","Reading ANI table\n"))
data<-read.table("ani_table",dec=".",sep="\t",row.names=NULL) #fastANI output table (all vs all)
colnames(data)<-c("Seq1","Seq2","ANI","Frags_used","Num_total_frags")
data2<-subset(data, ANI > 95) #Filter hits by 95% ANI, proposed threshold for viral species
cat(paste("[",Sys.time(),"]","Done\n"))

cat(paste("[",Sys.time(),"]","Reading length table\n"))
length<-read.table("Length_viruses.txt",header=F,sep="\t") #Viral sequence length
colnames(length)<-c("Virus","length")
cat(paste("[",Sys.time(),"]","Done\n"))

cat(paste("[",Sys.time(),"]","Merging tables\n"))
data2<-inner_join(data2, length, by = join_by(Seq1 == Virus), relationship = "many-to-many") #Merge fastANI and length tables by Seq1
data2<-inner_join(data2, length, by = join_by(Seq2 == Virus), relationship = "many-to-many") #Merge fastANI and length tables by Seq2
colnames(data2)<-c("Seq1","Seq2","ANI","Frags_used","Num_total_frags","Length_Seq1","Length_Seq2")
data2$Length_Seq1 <- as.numeric(data2$Length_Seq1)
data2$Length_Seq2 <- as.numeric(data2$Length_Seq2)
data2 <- data2 %>% mutate(Perc_genome_shared = ifelse((Length_Seq1 > Length_Seq2),(Frags_used*1000)/(Length_Seq2)*100,(Frags_used*1000)/(Length_Seq1)*100)) #Calculate % of genome shared based on Frags_used and genome length
data2$Perc_genome_shared[data2$Perc_genome_shared > 100] <- 100 #Due to genome fragmentation used by fastANI, sometimes the % genome shared is by artifact above 100. This line corrects this effect
cat(paste("[",Sys.time(),"]","Done\n"))


cat(paste("[",Sys.time(),"]","Building networks\n"))
edge_df <- data2[, c(1,2)] # Two column table with the connected IDs on each one


#Get names within a network
g2<-graph_from_data_frame(edge_df)
clusters1 = clusters(g2)
networks<-max(clusters1$membership) #Number of networks

cat(paste("[",Sys.time(),"]","Clusters identified:", networks, "\n"))
cat(paste("[",Sys.time(),"]","Get names per cluster\n"))

lista=list()
for (i in 1:networks){ 
  lista[[i]]<-V(g2)[clusters1$membership == i]
}
len<-length(lista) #Must be equal to variable "networks"

cat(paste("[",Sys.time(),"]","Identify each species\n"))

setDT(data2)

sp_identification <- function(i){
  df <- data2[Seq1 %in% names(lista[[i]]) | Seq2 %in% names(lista[[i]]), ]
  name<-paste("Species_", i,sep="")
  df$Species<-rep(name,each=nrow(df))
  return(df)
}

num_networks <- seq(1,length(lista), by=1)
res_sp <- lapply(num_networks, sp_identification)
df_final<-rbindlist(l=res_sp)

df_final_derep <- df_final %>% distinct() #Remove duplicated rows
df_final_derep2 <- subset(df_final_derep, Seq1 != Seq2) #Remove rows with self hits

min_num_genomes=20 # Minimum number of genomes per species. 20 for isolated viruses, 10 for uncultured viruses (fosmid and long-reads).
cat(paste("[",Sys.time(),"]","Identifying species with more than", min_num_genomes, "genomes","\n"))

setDT(df_final_derep2)
Num_strain_per_sp_min10<-data.frame(Species=character(0))
for (sp in unique(df_final_derep2$Species)){
  test <- df_final_derep2[Species == sp, ]
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

df_final_derep_min20 <- df_final_derep2 %>% filter(Species %in% Num_strain_per_sp_min10$Species)

cat(paste("[",Sys.time(),"]","Writing table\n"))
write.table(df_final_derep_min20, file="ani_table_final_sp_included.txt",sep="\t", col.names=F,row.names = F,quote = F)

cat(paste("[",Sys.time(),"]","All done\n"))
