library(gcplyr)
library(ggplot2)
library(heatwaveR)
library(pracma)
library(tidyverse)

# Function to process data
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


# Function to find peaks
peak_finder <- function(counts, mean_smooth){
  pre_peaks = findpeaks(counts$Smooth_counts) %>% as.data.frame()
  if (isempty(pre_peaks) == FALSE){
    index_peaks = pre_peaks %>% filter(V1 > mean_smooth)
    peaks = c()
    for (idx in index_peaks$V2){
      pk = counts$ANI_round[idx]
      peaks = append(peaks, pk)
    }
  } else {
    peaks <- 0
  }
  return(peaks)
}


# Function to find valleys
valley_finder <- function(counts, mean_smooth){
  pre_valleys = findpeaks(-counts$Smooth_counts) %>% as.data.frame()
  if (isempty(pre_valleys) == FALSE){
    index_valleys = pre_valleys %>% filter(-V1 < mean_smooth*0.5)
    valleys = c()
    for (idx in index_valleys$V2){
      vl = counts$ANI_round[idx]
      valleys = append(valleys, vl)
    }
  } else {
    valleys <- 0
  }
  return(valleys)
}


# Function to get range of valleys
get_range_valleys <- function(counts, mean_smooth, valleys){
  below_half_mean <- counts %>% filter(Smooth_counts < mean_smooth*0.5)
  
  #Get ranges of low frequency of pairs
  y=1
  low_freq_ranges <- list()
  for (x in seq(1:length(below_half_mean$ANI_round))){
    
    value = below_half_mean$ANI_round[x]
    
    if (x == 1){
      vec_values <- c()
      vec_values <- append(vec_values, value)
    }
    
    if (x > 1){
      test_eq = all.equal(value, (below_half_mean$ANI_round[x-1] + 0.1))
      if (test_eq == TRUE){
        vec_values <- append(vec_values, value)
      } else {
        low_freq_ranges[[y]] <- vec_values
        vec_values <- c()
        vec_values <- append(vec_values, value)
        y = y + 1
      }
    }
  }
  low_freq_ranges[[y]] <- vec_values
  
  # Check if valleys were found in area of low frequency of pairs
  accepted_valley_ranges <- list()
  for (z in 1:length(low_freq_ranges)){
    res_logical <- low_freq_ranges[z][[1]] %in% valleys
    out <- sum(res_logical, na.rm = TRUE)
    if (out > 0){
      accepted_valley_ranges <- append(accepted_valley_ranges, low_freq_ranges[z])
    }
  }
  # Return ranges of low frequency of pairs were valleys were also detected
  return(accepted_valley_ranges)
}


# Function to classify plots based on distribution
classify <- function(table_sp, peaks, range_valleys, area){
  avg_ANI <- mean(table_sp$ANI)
  
  # Check peaks
  check_peaks_in_area <- sum(peaks %in% area, na.rm=TRUE)
  if (check_peaks_in_area == 0){
    peaks_in_area <- FALSE
  } else {
    peaks_in_area <- TRUE
  }
  
  
  # Check if any range of valleys contains the values in the area of interest (99.2 - 99.8)
  index_pos_valleys <- list()
  for (m in seq(1:length(range_valleys))){
    if (sum(range_valleys[m][[1]] %in% area) > 0){
      index_pos_valleys <- append(index_pos_valleys, m)
    }
  }
  
  # Count the number of ranges of valleys with the values of interest
  if (length(index_pos_valleys) > 0){
    valley_in_area = TRUE
  } else {
    valley_in_area = FALSE
  }
  
  # Perform the classification
  if (avg_ANI < 99.8 & valley_in_area == TRUE){
    group = "Group 1"
  } else if (avg_ANI > 99.8){
    group = "Group 2"
  } else if (peaks_in_area == TRUE & valley_in_area == FALSE){
    group = "Group 4"
  } else {
    group = "Group 3"
  }
  
  return(list(group, index_pos_valleys, avg_ANI))
}


# Function to plot
plotting <- function(counts, peaks, valleys, sp, group){
  
  df_mean <- as.data.frame(mean(counts$Smooth_counts)*0.5)
  
  if (is.null(peaks) == FALSE){
    df_peaks <- as.data.frame(peaks)
  } else {
    df_peaks <- as.data.frame(0)
    colnames(df_peaks) <- c("peaks")
  }
  
  if (is.null(valleys) == FALSE){
    df_valleys <- as.data.frame(valleys)
  } else {
    df_valleys <- as.data.frame(0)
    colnames(df_valleys) <- c("valleys")
  }
  
  
  p <- ggplot(counts, aes(ANI_round, Num)) + geom_col(fill="lightblue", color="gray50") + 
    geom_flame(aes(y=mean_smooth*0.5, y2=Smooth_counts), alpha=0.4) +
    geom_line(aes(y=Smooth_counts, linetype="Smoothed", color="Smoothed"), show.legend = F) +
    geom_hline(data = df_mean, aes(yintercept = mean(counts$Smooth_counts)*0.5, linetype="Mean smoothed * 0.5", color="Mean smoothed * 0.5"), linewidth=1) +
    geom_vline(data = df_peaks, aes(xintercept = peaks, linetype="Peaks", color="Peaks"), linewidth=1.5, show.legend = F) + 
    geom_vline(data = df_valleys, aes(xintercept = valleys, linetype = "Valleys", color = "Valleys"), linewidth=1.5, show.legend = F) +
    scale_linetype_manual(name="Lines", values=c(1,3,2,4), guide = guide_legend(override.aes = list(color = c("black","black","blue","lightgreen"))), breaks = c("Smoothed", "Mean smoothed * 0.5", "Peaks", "Valleys")) +
    scale_color_manual(name="Lines", values = c("black","black","blue","lightgreen"), breaks = c("Smoothed", "Mean smoothed * 0.5", "Peaks", "Valleys")) +
    xlab("ANI (%)") + ylab("Frequency") + 
    labs(title = sp, subtitle = group) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0), limits = c(94.9,100.1)) +
    theme_bw() + theme(text=element_text(size=20), 
                       plot.title = element_text(size=18, colour = "gray20"),
                       plot.subtitle = element_text(size=15, colour = "gray30"),
                       legend.key.width = unit(2, "cm"))
  return(p)
}
  
# 1. Read
cat(paste("[",Sys.time(),"]","Reading table...\n"))
table <- read.table("ani_table_final_sp_included.txt", header=T, dec=".", sep="\t")
colnames(table) <- c("Seq1","Seq2","ANI","Frag_used","Frag_total","Length_Seq1","Length_Seq2","Perc_genome_shared","Species")

# 2. Process
df_final <- data.frame()
pdf("Species_classification_plots.pdf", width=16, height=9)

sp_unique <- unique(table$Species)

cat(paste("[",Sys.time(),"] Looping over ", length(sp_unique), " species...\n", sep = ""))
for (sp in sp_unique){
  table_sp <- table %>% filter(Species == sp)
  counts = process(table_sp)
  
  # 3. Find peaks and valleys
  if ('Smooth_counts' %in% names(counts)){
    mean_smooth = mean(counts$Smooth_counts)
  } else {
    cat(paste("Smoothed data absent. Skipping species", sp, "\n",sep=" "))
    next
  }

  peaks = peak_finder(counts, mean_smooth)
  
  valleys = valley_finder(counts, mean_smooth)
  
  # 4. Find range of valleys
  range_valleys <- get_range_valleys(counts, mean_smooth, valleys)
  
  # 5. Classify distribution
  area <- round(seq(99.2, 99.8, by=0.1),1) # The gap is the area between 99.2 - 99.8% ANI
  list_res <- classify(table_sp, peaks, range_valleys, area)
  group <- list_res[[1]]
  index_pos_valleys <- list_res[[2]]
  avg_ANI <- list_res[[3]]
  
  # 6. Plot
  p <- plotting(counts, peaks, valleys, sp, group)
  print(p)
  
  # 7. Print result
  if (length(index_pos_valleys) > 0){
    valleys_to_print <- toString(paste0(range_valleys[[index_pos_valleys[[1]]]]))
  } else {
    valleys_to_print <- "-"
  }
  
  final_res <- paste(sp, avg_ANI, toString(paste0(peaks)), valleys_to_print, group, "\n",sep="\t")
  cat(final_res)
  df_final <- rbind(df_final, c(sp, avg_ANI, toString(paste0(peaks)), valleys_to_print, group))
}

# 8. Close pdf
dev.off()

# 9. Provide column names and save output to a tab-delimited file
colnames(df_final) <- c("Species", "Avg_ANI", "Peaks", "Valleys_range", "Group")
write.table(df_final, "Peak_valleys_per_sp.txt", col.names = T, row.names = F, sep="\t", quote=F)

