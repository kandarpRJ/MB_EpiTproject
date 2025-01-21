library(rtracklayer)
library(dplyr)

setwd("~/data/RNA_modifications/MB_epi")

base_path <- "/media/pereralab/SSD/MB_epitranscriptomics/pileups_chm13v2"

sample_names <- list.files(path = base_path, pattern = ".a.cov4_pmod10_mod2_omod0.bed")
sample_names <- gsub(".pileup.a.cov4_pmod10_mod2_omod0.bed", "", sample_names)
mod_names <- c("a", "p", "m", "i")

count_modifications <- function(sample_names, mod_names, base_path) {
  count_table <- data.frame(Sample = character(),
                            Modification = character(),
                            Count = integer(),
                            stringsAsFactors = FALSE)

  for (sample in sample_names) {
    for (mod in mod_names) {
      bed_file <- file.path(base_path, paste0(sample, ".pileup.", mod, ".cov4_pmod10_mod2_omod0.bed"))
      
      if (file.exists(bed_file)) {
        mod_sites <- read.table(bed_file, sep = "\t", row.names = NULL, header = F)
        mod_sites <- mod_sites[,1:3]
        colnames(mod_sites)<-c("chrom", "start", "end")
        mod_sites <- as(mod_sites, "GRanges")
        mod_sites <- GRanges(mod_sites)
        
        # Count number of modifications
        mod_count <- length(mod_sites)
      } else {
        mod_count <- 0  # If file is missing, set count to 0
      }
      
      # Store results in the dataframe
      count_table <- rbind(count_table, data.frame(Sample = sample, 
                                                   Modification = mod, 
                                                   Count = mod_count))
    }
  }
  
  return(count_table)
}

modification_counts <- count_modifications(sample_names, mod_names, base_path)

modification_pivot <- modification_counts %>%
  pivot_wider(names_from = Modification, values_from = Count, values_fill = 0)

print(modification_pivot)




