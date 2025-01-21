library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/data/RNA_modifications/MB_epi")
# Function to count transcript biotypes from GFF3 and BED files
calculate_biotype_counts <- function(gff3_file, bed_files, mod_names) {
  
  # Read the GFF3 file
  gff_data <- import.gff3(gff3_file)
  
  # Filter for transcript features and extract biotypes
  transcripts <- gff_data[gff_data$type == "transcript"]
  transcript_biotypes <- transcripts$transcript_biotype
  
  # Create an empty list to store results
  result_list <- list()
  
  # Loop through each BED file and count biotypes
  for (i in seq_along(bed_files)) {
    bed_file <- bed_files[i]
    mod_name <- mod_names[i]
    
    # Import BED file containing RNA modifications
    mod_sites <- import.bed(bed_file)
    
    # Find overlaps between modifications and transcripts
    overlaps <- findOverlaps(transcripts, mod_sites)
    
    # Get overlapping transcript biotypes
    overlapping_biotypes <- transcript_biotypes[unique(queryHits(overlaps))]
    
    # Count biotypes
    biotype_counts <- table(overlapping_biotypes)
    biotype_df <- as.data.frame(biotype_counts)
    colnames(biotype_df) <- c("biotype", "count")
    biotype_df$modification <- mod_name
    
    result_list[[mod_name]] <- biotype_df
  }
  
  # Combine all data into one dataframe
  combined_df <- bind_rows(result_list)
  
  return(combined_df)
}

# Define input files
gff3_file <- "~/ref/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"
bed_files <- c("/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_m6A.cov4_pmod10_mod2_omod0.bed",
               "/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_m5C.cov4_pmod10_mod2_omod0.bed",
               "/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_I.cov4_pmod10_mod2_omod0.bed",
               "/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2/all_psU.cov4_pmod10_mod2_omod0.bed")
mod_names <- c("m6A", "m5C", "Inosine", "pseudoU")

# Run the function to get biotype counts
biotype_counts_df <- calculate_biotype_counts(gff3_file, bed_files, mod_names)
biotype_counts_df$count<-log(biotype_counts_df$count)

# Plot combined bar graph
ggplot(biotype_counts_df, aes(x = biotype, y = count, fill = modification)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Transcript Biotype Distribution for RNA Modifications",
       x = "Transcript Biotype",
       y = "log(Count)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18, face = "bold"),  
        axis.title.y = element_text(size = 18, face = "bold"),  
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +  
  scale_fill_manual(values = c("m6A" = "darkred", 
                               "m5C" = "darkblue", 
                               "Inosine" = "darkgreen", 
                               "pseudoU" = "darkorange")) 

# Save the plot
ggplot2::ggsave("all_RNA_mod_biotype_distribution.png", width = 13, height = 6)
