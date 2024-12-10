args <- commandArgs(trailingOnly=TRUE)

library(MOSAIC)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(viridis)

#Directory containing MOSAIC output and snp files
mosaic_input_dir = args[1]
#Sample name
SAMPLE= args[2]
#Directory to save output files
out_dir= args[3]
#Path to data
mosaic_data= args[4]

#Global variables
Anc0="European"
Anc1="PreContactDingo"
threshold=0.9

la_results = paste0(mosaic_input_dir,"localanc_",mosaic_data)
model_results = paste0(mosaic_input_dir,mosaic_data)

#load files
load(model_results)
load(la_results)

# localanc gives the local ancestry at each grid point
# get local ancestry probabilities at each SNP
local_pos=grid_to_pos(localanc, mosaic_input_dir, g.loc, chrnos)

for(chromosome in 1:38) {
  local_pos_df <- reshape2::melt(local_pos[[chromosome]])
  colnames(local_pos_df) <- c("ancestry", "haplotype", "number", "prob")
  
  snp=read.table(paste0(mosaic_input_dir,"snpfile.",chromosome))
  snp_row=nrow(snp)
  sequence_length <- snp[snp_row,4]
  
  #create dummy list
  datalist_mosaic_hapl1 = list()
  datalist_mosaic_hapl2 = list()
  
  for(i in seq(1, max(local_pos_df$haplotype), 2)){
    
    mosaic_sample_hap1 <- local_pos_df %>%
      filter(haplotype %in% c(i)) %>% 
      pivot_wider(names_from = ancestry, values_from = prob, names_prefix="prob_anc") %>%
      mutate(AN1 = ifelse(prob_anc1>prob_anc2, Anc0, Anc1)) %>%
      mutate(AN1 = ifelse(prob_anc1<threshold & prob_anc2<threshold, 'pop_no', AN1)) %>%
      mutate(sample = paste(SAMPLE, (i+1)/2, sep = ""))
    mosaic_sample_hap1_snp <- cbind(mosaic_sample_hap1, snp[,4, drop=FALSE])
    
    datalist_mosaic_hapl1[[i]] <- mosaic_sample_hap1_snp
    
    mosaic_sample_hap2 <- local_pos_df %>%
      filter(haplotype %in% c(i+1)) %>% 
      pivot_wider(names_from = ancestry, values_from = prob, names_prefix="prob_anc") %>%
      mutate(AN2 = ifelse(prob_anc1>prob_anc2, Anc0, Anc1)) %>%
      mutate(AN2 = ifelse(prob_anc1<threshold & prob_anc2<threshold, 'pop_no', AN2)) %>%
      mutate(sample = paste(SAMPLE, (i+1)/2, sep = ""))
    mosaic_sample_hap2_snp <- cbind(mosaic_sample_hap2, snp[,4, drop=FALSE])
    
    datalist_mosaic_hapl2[[i]] <- mosaic_sample_hap2_snp
  }
  
  big_data_mosaic_hap1 = do.call(rbind, datalist_mosaic_hapl1) #combine all together hapl1
  big_data_mosaic_hap2 = do.call(rbind, datalist_mosaic_hapl2) #combine all together hapl2
  
  
  #create duplicated first row and replace with site=1
  big_data_mosaic_hap1_dupl <- big_data_mosaic_hap1 %>%
    group_by(sample) %>%
    dplyr::slice(1L, row_number()) %>% #replicate the first row
    mutate(V4 = ifelse(row_number() == 1, 0, V4)) #replace the first snp to be 0
  
  big_data_mosaic_hap2_dupl <- big_data_mosaic_hap2 %>%
    group_by(sample) %>%
    dplyr::slice(1L, row_number()) %>% #replicate the first row
    mutate(V4 = ifelse(row_number() == 1, 0, V4)) #replace the first snp to be 0
  
  # Merge datasets based on the V4 column
  merged_df <- merge(big_data_mosaic_hap1_dupl, big_data_mosaic_hap2_dupl, by = "V4")
  
  filtered_df <- merged_df %>%
    select(-number.x, -prob_anc2.x, -haplotype.y, -number.y, -prob_anc2.y, -sample.y) 
  colnames(filtered_df) <- c("SNP", "chr", "prob_European_1", "ANC1", "Sample", "prob_European_2", "ANC2")
  
  final_df <- filtered_df %>% select(Sample, chr, SNP, prob_European_1, ANC1, prob_European_2, ANC2)
  
  write.table(final_df, file=paste0(out_dir,SAMPLE,"_",chromosome,"_anc"), quote=FALSE, sep='\t', row.names=FALSE)
}  
