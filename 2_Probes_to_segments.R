# 
# Probes to segments
#

# Define gain/loss for each segment and gather the segments of all samples in one data frame

# Load folder names ------------------------------------------------------------
library(data.table)
library(DescTools)
library(ggplot2)
library(tidyr)
library(dplyr)
library(feather)
library(tidyverse)

setwd("D:/SNP 6/Rohdaten/cel")
folders <- grep(list.files(), pattern = '\\.CEL|\\.png|\\.ARR|\\.pdf|\\.csv|\\CHP|preprocessed', 
                inv = TRUE, value = TRUE)

# Loop for probe to segments ---------------------------------------------------

probe_zuweiser <- function(segments, probes){
  # get probe subset with correct chromosome
  probes.tmp <- probes[probes$Chromosome == segments[1], ]
  # get probes in segment
  probes.val <- probes.tmp[as.numeric(probes.tmp$Start) >= as.numeric(segments[2]) & 
                             as.numeric(probes.tmp$End) <= as.numeric(segments[3]), ]
  if (nrow(probes.val) != 0){
    probes.val$LogR <- as.numeric(segments[5])
    probes.val$AimBa <- as.numeric(segments[7])
  } else probes.val <- c(NA, NA)
  return(probes.val)
}
for( i in 1:length(folders)){
  # Load Probes data
  load(paste0(folders[i], "/rawcopy.Rdata"))
  
  # get probes Id and position
  probes <- probes.txt[, c("Probe.Set.ID", "Chromosome", "Start", "End", 
                           "Value", "smoothed", "B.Allele.Frequency")]
  names(probes)[5:7] <- c("Probes.txt.value", "Probes.txt.smoothed", "Probes.txt.BAF")                         
  probes <- probes[probes$Chromosome != "chrM",]
  probes$LogR <- NA
  probes$AimBa <- NA
  
  seg <- fread(paste0(folders[i],"/segments2.txt"))
  seg <- as.data.frame(seg[, -11])
  seg <- seg[seg$Chromosome != "chrM",]
  
  
  
  values <- apply(seg, 1, probe_zuweiser, probes = probes)
  values <- do.call(rbind, values)
  values$status <- "default"
  values$status[values$LogR > 0.2] <- "gain"
  values$status[values$LogR < -0.1] <- "loss"
  saveRDS(values, paste0(folders[i], "/probe.values.rds"))
  saveRDS(table(values[, c("Chromosome", "status")]),
          paste0(folders[i], "/gl_table.rds"))
  
}

# Load all Segment values in a data frame --------------------------------------
pheno <-  read.csv("2018 10 10 Wien PRCa Update_prepared.csv", 
                   stringsAsFactors = FALSE)
folder.files <- sapply(folders, 
                       function(x){
                         x <- strsplit(x, "_", fixed = TRUE)[[1]][3]
                         return(x)
                       })
folder.files[208:length(folders)] <- sapply(folders[208:length(folders)],
                                            function(x){
                                              x <- strsplit(x, "_", fixed = TRUE)[[1]][2]
                                              return(x)
                                            })

seg <- fread(paste0(folders[1],"/segments2.txt"))
seg <- as.data.frame(seg[, -11])
seg <- seg[seg$Chromosome != "chrM",]
seg$Proband <- names(seg)[5]
names(seg)[5] <- "LogR"
names(seg)[7] <- "Allelic.Imbalance"
seg$status <- "default"
seg$status[seg$LogR > 0.6] <- "gain"
seg$status[seg$LogR < -0.6] <- "loss"
seg$cc <- pheno[pheno$File_short == folder.files[1], "status"]


segments_df <- seg
for( i in 2:length(folders)){
  cat(i, "of 457 processed\n")
  seg <- fread(paste0(folders[i],"/segments2.txt"))
  seg <- as.data.frame(seg[, -11])
  seg <- seg[seg$Chromosome != "chrM",]
  seg$Proband <- names(seg)[5]
  names(seg)[5] <- "LogR"
  names(seg)[7] <- "Allelic.Imbalance"
  seg$status <- "default"
  seg$status[seg$LogR > 0.6] <- "gain"
  seg$status[seg$LogR < -0.6] <- "loss"
  if(sum(pheno$File_short == folder.files[i]) == 0) next
  seg$cc <- pheno[pheno$File_short == folder.files[i], "status"]
  segments_df <- rbind(segments_df, seg)
}
feather::write_feather(segments_df,
                       "all_probands_segments_df.feather")
