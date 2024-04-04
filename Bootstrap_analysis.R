library(ggplot2)
library(pracma)
library(gcplyr)
library(tidyverse)
library(data.table)

# On Mac
setwd("/Users/briquelme3/Dropbox (GaTech)/ANI_virus/fastANI/Prokaryotic_phages/IMG_VR")

# Function to subsample the same number of pairs per species
subsample <- function(table, sp_unique){
  num_pairs = 150
  df_subsampled <- data.frame()
  for (sp in sp_unique){
    table_sp <- table[Species == sp]
    if (length(table_sp$ANI) > num_pairs){
      random_data<-sample_n(table_sp, num_pairs)
      df_subsampled <- rbind(df_subsampled, random_data)
    }
  }
  return(df_subsampled)
}


# Function to find peaks
peak_finder <- function(counts, mean_smooth){
  pre_peaks = findpeaks(counts$Smooth_counts) %>% as.data.frame()
  if (isempty(pre_peaks) == FALSE){
    index_peaks = pre_peaks %>% filter(V1 > mean_smooth)
    peaks = c()
    for (idx in 1:length(index_peaks$V1)){
      index_peaks$ANI[idx] <- counts$ANI_round[index_peaks$V2[idx]] 
    }
  } else {
    peaks <- 0
  }
  return(index_peaks)
}


# Function to find valleys
valley_finder <- function(counts, mean_smooth){
  pre_valleys = findpeaks(-counts$Smooth_counts) %>% as.data.frame()
  if (isempty(pre_valleys) == FALSE){
    index_valleys = pre_valleys %>% filter(-V1 < mean_smooth*0.5)
    index_valleys$V1 <- -index_valleys$V1
    for (idx in 1:length(index_valleys$V1)){
      index_valleys$ANI[idx] <- counts$ANI_round[index_valleys$V2[idx]] 
    }
  } else {
    index_valleys <- 0
  }
  return(index_valleys)
}


# Function to process
process <- function(table_sp){
  table_sp$ANI_round <- round(table_sp$ANI,1)
  counts <- table_sp %>% group_by(ANI_round) %>% summarise(Num=n())
  
  min_bin <- min(counts$ANI_round)
  
  for (num in seq(min_bin,100,by=0.1)){
    if (round(num,1) %in% counts$ANI_round == FALSE){
      counts <- rbind(counts,c(num, 0))
    }
  }
  
  counts <- counts[order(counts$ANI_round),]
  
  tryCatch(
    expr = {
      counts$Smooth_counts <- smooth_data(y = counts$Num, x = counts$ANI_round, sm_method = "gam")
      counts$cumulative <- cumsum(counts$Num)
    },
    error = function(e){
      cat("Smooth step raised an error. ")
    })
  
  return(counts)
}

cat(paste("[",Sys.time(),"]","Reading table...\n"))
table <- fread("ani_table_final_sp_included.txt", header=T, dec=".", sep="\t")
colnames(table) <- c("Seq1","Seq2","ANI","Frag_used","Frag_total","Length_Seq1","Length_Seq2","Perc_genome_shared","Species")

df_final_valley_list <- data.frame()
df_final_peak_list <- data.frame()

# Bootstrap subsampling
x=0
while (x < 1000){
  sp_unique <- unique(table$Species)
  cat(paste("[",Sys.time(),"] Subsampling ", length(sp_unique), " species...\n", sep = ""))
  df_subsampled <- subsample(table, sp_unique)
  
  cat(paste("[",Sys.time(),"] Smoothing data...\n", sep = ""))
  counts <- process(df_subsampled)
  
  cat(paste("[",Sys.time(),"] Finding peaks and valleys...\n", sep = ""))
  mean_smooth = mean(counts$Smooth_counts)
  peaks <- peak_finder(counts, mean_smooth)
  valleys <- valley_finder(counts, mean_smooth)
  
  df_final_peak_list <- rbind(df_final_peak_list, peaks)
  df_final_valley_list <- rbind(df_final_valley_list, valleys)
  
  x=x+1
  print(x)
}

# Summarize peaks and valleys identified in the bootstrap subsampling
summary_peak <- df_final_peak_list %>% 
  group_by(ANI) %>% summarise(count=n(), Height=mean(V1)) %>% add_column(Type="Peak")
colnames(summary_peak) <- c("ANI", "Count", "Height", "Type")

summary_valley <- df_final_valley_list %>%
  group_by(ANI) %>% summarise(count=n(), Height=mean(V1)) %>% add_column(Type="Valley")
colnames(summary_valley) <- c("ANI","Count","Height","Type")

merged_peaks_valleys <- rbind(summary_valley, summary_peak)

# Plot smoothed counts
counts <- counts %>% mutate(Type=ifelse(Smooth_counts < mean_smooth*0.5, "Valley", "Peak"))

a <- ggplot(counts, aes(x=ANI_round, y=Smooth_counts, fill=Type, color=Type)) + #geom_col(fill="lightblue", color="darkblue") +
  geom_col() + scale_fill_manual(values = c("lightblue", "#31e47d")) + scale_color_manual(values = c("darkblue", "#01863a")) +
  xlab("ANI (%)") + ylab("Smoothed\ncount") + 
  scale_x_continuous(limits=c(94.9,100.1)) +
  geom_hline(aes(yintercept = mean_smooth * 0.5, linetype = "Smoothed mean * 0.5"), linewidth=1) +
  scale_linetype_manual(values = c(3)) +
  theme_bw() + theme(text=element_text(size=25), axis.text.x = element_blank(), 
                     axis.title.x = element_blank(), legend.position = c(0.15, 0.75),
                     legend.title = element_blank(), legend.key.size = unit(1,"cm"))

# Plots peaks and valleys find on each bootstrap
b <- ggplot(merged_peaks_valleys, aes(x=ANI,y=Count, color=Type, fill=Type)) + geom_col() + 
  xlab("ANI (%)") + ylab("Number\nbootstraps") +
  scale_fill_manual(values = c("lightblue", "#31e47d")) + scale_color_manual(values = c("darkblue", "#01863a")) +
  scale_x_continuous(limits=c(94.9,100.1)) +
  theme_bw() + theme(text=element_text(size=25), legend.position = "None")
b

# Merge plots
pdf("Bootstrap_plot.pdf", width=13, height=13)
print(ggpubr::ggarrange(a, b, ncol=1))
dev.off()

