library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggpubr)

plotting <- function (tabla, species) {
  tabla['ANI_round'] <- round(tabla$ANI,2)
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
  
  a<-ggplot(ani_round_freq, aes(x=ANI_round,y=n)) + geom_col(color="black",fill="#cba4d8") + 
    scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0), limits=c(94.9,100.1)) +
    xlab("ANI (%)") + ylab("Frequency") +
    theme_bw() + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  
  perc_round_freq<-tabla %>% group_by(Perc_round) %>% summarise(n=n())
  
  l<-ggplot(perc_round_freq, aes(x=Perc_round,y=n)) + geom_col(color="black",fill="#cba4d8",width=0.9) + 
    scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0), limits=c(round(min(tabla$Perc_genome_shared),0)-1.1,100.5)) +
    xlab("% shared genome") + ylab("Frequency") + 
    theme_bw() + theme(axis.text.y  = element_blank(), axis.title.y = element_blank())
  
  pdf(paste("Plot_",species,".pdf",sep=""), width=14,height=9)
  print(ggarrange(a,NULL, m, l + coord_flip(), ncol = 2, nrow = 2, heights = c(0.25,1), widths  = c(1,0.3)))
  dev.off()
}


#setwd() if needed

data<-read.table("ani_table_CDC_USA_final",dec=".",sep="\t",row.names=NULL) #ANI table obtained with FastANI
colnames(data)<-c("Seq1","Seq2","ANI","Frags_used","Num_total_frags")
data2<-subset(data, ANI > 95) #Remove pairs sharing less than 95% ANI

length<-read.table("Length_1250_genomes_each_VOC_variant_CDC_USA.txt",header=F,sep="\t")
colnames(length)<-c("Virus","length")
data2<-inner_join(data2, length, by = join_by(Seq1 == Virus))
data2<-inner_join(data2, length, by = join_by(Seq2 == Virus))
colnames(data2)<-c("Seq1","Seq2","ANI","Frags_used","Num_total_frags","Length_Seq1","Length_Seq2")
data2["Species"] <- c("Sars-Cov2")
data2 <- data2 %>% mutate(Perc_genome_shared = ifelse((as.integer(Length_Seq1) > as.integer(Length_Seq2)),(as.integer(Frags_used)*1000)/(as.integer(Length_Seq2))*100,(as.integer(Frags_used)*1000)/(as.integer(Length_Seq1))*100))
data2$Perc_genome_shared[data2$Perc_genome_shared > 100] <- 100

df_final_derep <- data2 %>% distinct() #Remove duplicated rows
df_final_derep2 <- subset(df_final_derep, Seq1 != Seq2) #Remove rows with self hits
df_final_derep2 <- subset(df_final_derep2, Num_total_frags > 20) #Remove rows with viruses shorter than 20 kb

#Plot ANI vs Percentage genome shared
npairs=1000000
sp="SARS-CoV-2"
plot_data <- df_final_derep2 %>% slice_sample(n=npairs)
plotting(plot_data, sp)
  

#WHO variant information and histograms
metadata <- read.table("Metadata_1250_genomes_each_VOC_variant_CDC_USA.txt",header=F,sep="\t")
metadata_variant <- metadata[, c('V1', 'V6') ]
colnames(metadata_variant) <- c("Accession","Variant")

df_final_derep3 <- inner_join(df_final_derep2, metadata_variant ,by = join_by(Seq1 == Accession))
df_final_derep3 <- inner_join(df_final_derep3, metadata_variant ,by = join_by(Seq2 == Accession))
colnames(df_final_derep3)<-c("Seq1","Seq2","ANI","Frags_used","Num_total_frags","Length_Seq1","Length_Seq2", "Species","Perc_genome_shared","Variant_Seq1","Variant_Seq2")

df_final_derep3 <- df_final_derep3 %>% mutate(Variant_comparison=ifelse(Variant_Seq1 == Variant_Seq2, "Same variant","Different variants")) 

a<-ggplot(df_final_derep3, aes(x=ANI)) + geom_histogram(color="black",fill="#ebdbb0", alpha=0.6, bins=100) + 
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0), limits=c(99.4,100.05), breaks=seq(98,100,by=0.1)) +
  xlab("ANI (%)") + ylab("Frequency") +
  geom_vline(xintercept = 99.83, linetype="dashed", linewidth=1, colour="gray50") +
  theme_bw() + theme(text=element_text(size=20))

b<-ggplot(df_final_derep3, aes(x=ANI, fill=Variant_comparison)) + 
  geom_histogram(data=subset(df_final_derep3, Variant_comparison == "Same variant"),aes(color="Same variant"),bins=100, alpha=0.5, color="black") +
  geom_histogram(data=subset(df_final_derep3, Variant_comparison == "Different variants"),aes(color="Different variants"),bins=100, alpha=0.5, color="black") +
  scale_fill_manual(name="Variant comparison", values = c("Same variant" = "#619CFF", "Different variants" = "#F8766D")) +
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0), limits=c(99.4,100.05), breaks=seq(99,100,by=0.1)) +
  geom_vline(xintercept = 99.83, linetype="dashed", linewidth=1, colour="gray50") +
  xlab("ANI (%)") + ylab("Frequency") +
  theme_bw() + theme(text=element_text(size=20), legend.position = c(0.13,0.82))

ggarrange(a,b,ncol=1,nrow=2, labels = c("A","B"))




#Numeric data
threshold <- 99.84 # ANI bin with the lowest number of pairs between the two modes (peaks)
##Above threshold
res <- df_final_derep3 %>% filter(ANI > threshold) %>% 
  group_by(Variant_comparison) %>% summarise(num=n()) %>% mutate(Perc= round(num / sum(num) * 100, 1)) %>%
  pull(Perc)
cat(paste(res[1],"% of comparisons ABOVE ", threshold, "% ANI occur between different variants and ",res[2], "% between the same variant","\n",sep="")) #res[1] es el % de comparaciones entre distintas variantes, res[2] es el % de comparaciones entre la misma variante

##Below threshold
res <- df_final_derep3 %>% filter(ANI < threshold) %>% 
  group_by(Variant_comparison) %>% summarise(num=n()) %>% mutate(Perc= round(num / sum(num) * 100, 1)) %>%
  pull(Perc)
cat(paste(res[1],"% of comparisons BELOW ", threshold, "% ANI occur between different variants and ",res[2], "% between the same variant","\n",sep="")) #res[1] es el % de comparaciones entre distintas variantes, res[2] es el % de comparaciones entre la misma variante

