#' @title miRNA target prediction by position of prediction
#' @description Combine results from RNAhybrid and miRanda by position, remove target isoforms with identical 3'UTR seq, and filter for internal/bulge loop, return summary of results
#' A bug was found (on 20210823) for results of RNAhybrid: for every 1000 bp of the 3'UTR,
#' the predicted miRNA binding sites location will be decreased by 1 bp
#' No such bug was found in Dme
#' This script attempted to correct for the bug, then follow standard data analysis ("miRNA_target_prediction_pipeline.R")


org = "Tco" # c("Dme", "Tco", "Hho", "Nda")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(here::here())
setwd(paste0("../", org))

# load required library
library(xlsx)
library(tidyverse)
library(Biostrings)
library(stringr)
library(dplyr)

{
  # Read result files of RNAhybrid and miRanda (.xls)
  RNAhybrid_217 <- read.delim("Tco_Cody.3utr.fa.RNAhybrid_targets_seed-region_2-7_12-17.xls", stringsAsFactors=FALSE, sep="\t")
  RNAhybrid_217 <- read.delim("Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.RNAhybrid_targets_seed-region_2-7_12-17.xls", stringsAsFactors=FALSE, sep="\t")
  RNAhybrid_217 <- read.delim("Dme_3_prime_UTR_dedup.fasta.RNAhybrid_targets_seed-region_2-7_12-17.xls", stringsAsFactors=FALSE, sep="\t")
  head(RNAhybrid_217)
  miranda <- read.delim("Tco_Cody.3utr.fa.miranda_with_strict.hit.xls", stringsAsFactors=FALSE, sep="\t")
  miranda <- read.delim("Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.miranda_with_strict.hit.xls", stringsAsFactors=FALSE, sep="\t")
  miranda <- read.delim("Dme_3_prime_UTR_dedup.fasta.miranda_with_strict.hit.xls", stringsAsFactors=FALSE, sep="\t")
  head(miranda)
}

validate_length <- function(RNAhybrid_217, miranda) {
  RNAhybrid_target_unique <- RNAhybrid_217[!duplicated(RNAhybrid_217$target_id),]
  miranda_target_unique <- miranda[!duplicated(miranda$target),]
  
  RNAhybrid_miranda_target_unique <- merge(RNAhybrid_target_unique, miranda_target_unique, by.x="target_id", by.y="target")
  head(RNAhybrid_miranda_target_unique)
  print(RNAhybrid_miranda_target_unique$target_length == RNAhybrid_miranda_target_unique$target_Len)
  # not all the 3' UTR length are the same after merging the two results
  RNAhybrid_miranda_target_unique <- 
    RNAhybrid_miranda_target_unique[!RNAhybrid_miranda_target_unique$target_length == RNAhybrid_miranda_target_unique$target_Len, ]
  return(RNAhybrid_miranda_target_unique)
}
RNAhybrid_miranda_target_unique <- validate_length(RNAhybrid_217, miranda)

# Positions correction for RNAhybrid
RNAhybrid_correct <- function(RNAhybrid_217) {
  RNAhybrid_217_positioned <- RNAhybrid_217
  # RNAhybrid_217_positioned <- RNAhybrid_217[-which(RNAhybrid_217$miRNA_aln == ""),]
  # length correction
  RNAhybrid_217_positioned$target_length_seed <- 
    ifelse(nchar(RNAhybrid_217_positioned$target_length)<=3, 0, 
    as.numeric(substr(RNAhybrid_217$target_length, 1, nchar(RNAhybrid_217_positioned$target_length)-3)))
  
  RNAhybrid_217_positioned$target_length <- 
    as.numeric(RNAhybrid_217_positioned$target_length)
  
  RNAhybrid_217_positioned$target_length_corrected <- 
    RNAhybrid_217_positioned$target_length + RNAhybrid_217_positioned$target_length_seed
  
  RNAhybrid_217_positioned$target_length_corrected_seed <- 
    ifelse(nchar(RNAhybrid_217_positioned$target_length_corrected)<=3, 0, 
    as.numeric(substr(RNAhybrid_217_positioned$target_length_corrected, 1, nchar(RNAhybrid_217_positioned$target_length_corrected)-3)))
  
  RNAhybrid_217_positioned$target_length_seed_match <-
    RNAhybrid_217_positioned$target_length_seed == RNAhybrid_217_positioned$target_length_corrected_seed
  
  RNAhybrid_217_positioned <- 
    RNAhybrid_217_positioned %>% 
    relocate(c("target_length_seed", "target_length_corrected", "target_length_corrected_seed", "target_length_seed_match"), .after=target_length)
  
  # predicted sites position correction
  RNAhybrid_217_positioned$position_seed <- 
    ifelse(nchar(RNAhybrid_217_positioned$position)<=3, 0, 
    as.numeric(substr(RNAhybrid_217_positioned$position,1,nchar(RNAhybrid_217_positioned$position)-3)))
  
  RNAhybrid_217_positioned$position_corrected <- 
    RNAhybrid_217_positioned$position + RNAhybrid_217_positioned$position_seed
  
  RNAhybrid_217_positioned$position_corrected_seed <- 
    ifelse(nchar(RNAhybrid_217_positioned$position_corrected)<=3, 0, 
    as.numeric(substr(RNAhybrid_217_positioned$position_corrected,1,nchar(RNAhybrid_217_positioned$position_corrected)-3)))
  
  RNAhybrid_217_positioned$position_corrected_seed_match <-
    RNAhybrid_217_positioned$position_seed == RNAhybrid_217_positioned$position_corrected_seed
  
  RNAhybrid_217_positioned <- 
    RNAhybrid_217_positioned %>% relocate(c("position_seed","position_corrected","position_corrected_seed","position_corrected_seed_match"),.after=position)
  
  return(RNAhybrid_217_positioned)
}
RNAhybrid_217_positioned <- RNAhybrid_correct(RNAhybrid_217)

# RNAhybrid: Concatenate "target_id" and "miRNA_id" and "target length" and "Positions"
RNAhybrid_concat <- function(RNAhybrid_217_positioned) {
  miRNA_structure <- paste(RNAhybrid_217_positioned$target_aln, RNAhybrid_217_positioned$target_hit, 
                           RNAhybrid_217_positioned$miRNA_hit, RNAhybrid_217_positioned$miRNA_aln,sep="\n")
  RNAhybrid_217_positioned_structured <- cbind(RNAhybrid_217_positioned,miRNA_structure)
  
  RNAhybrid_217_positioned_structured$target_mirna <- 
    paste(RNAhybrid_217_positioned_structured$target_id, RNAhybrid_217_positioned_structured$miRNA_id,sep="_")
  RNAhybrid_217_positioned_structured$target_mirna_len <- 
    paste(RNAhybrid_217_positioned_structured$target_mirna, paste0("len",RNAhybrid_217_positioned_structured$target_length_corrected), sep="_")
  RNAhybrid_217_positioned_structured$target_mirna_pos<- 
    paste(RNAhybrid_217_positioned_structured$target_mirna, paste0("pos",RNAhybrid_217_positioned_structured$position_corrected), sep="_")
  RNAhybrid_217_positioned_structured$target_mirna_len_pos <- 
    paste(RNAhybrid_217_positioned_structured$target_mirna_len, paste0("pos",RNAhybrid_217_positioned_structured$position_corrected), sep="_")
  
  return(RNAhybrid_217_positioned_structured)
}
RNAhybrid_217_positioned_structured <- RNAhybrid_concat(RNAhybrid_217_positioned)

write.csv(RNAhybrid_217_positioned_structured, "Tco_Cody.3utr.fa.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", row.names=FALSE)
write.csv(RNAhybrid_217_positioned_structured, "Hho_Cody.3utr.fa.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", row.names=FALSE)
write.csv(RNAhybrid_217_positioned_structured, "Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", row.names=FALSE)
write.csv(RNAhybrid_217_positioned_structured, "Dme_3_prime_UTR_dedup.fasta.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", row.names=FALSE)
####

# Read corrected files of RNAhybrid 
RNAhybrid_217_positioned_structured <- read.csv("Tco_Cody.3utr.fa.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", fileEncoding="UTF-8-BOM")
RNAhybrid_217_positioned_structured <- read.csv("Hho_Cody.3utr.fa.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", fileEncoding="UTF-8-BOM")
RNAhybrid_217_positioned_structured <- read.csv("Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", fileEncoding="UTF-8-BOM")
RNAhybrid_217_positioned_structured <- read.csv("Dme_3_prime_UTR_dedup.fasta.RNAhybrid_targets_seed-region_2-7_12-17_positioned_structured.xls.csv", fileEncoding="UTF-8-BOM")

head(RNAhybrid_217_positioned_structured)


# Positions expansion for miRanda (divide multiple "Positions" into separate rows and data frames)
head(miranda)
{
  miranda$Positions <- str_squish(miranda$Positions)
  miranda$Positions_Gap <- str_count(miranda$Positions," ")
  
  miranda$target_mirna <- paste(miranda$target, sub("^>>", "", miranda$miRNAs), sep="_") #formatting
  miranda$target_mirna_len <- paste(miranda$target_mirna, paste0("len", miranda$target_Len), sep="_")

  # Divide into two data frames based on "Positions_Gap"
  miranda_pos_1 <- miranda[miranda$Positions_Gap == 0, ]
  miranda_pos_not1 <- miranda[miranda$Positions_Gap > 0, ]
}

# miRanda: generate column "target_mirna_pos" / "target_mirna_len_pos"
{
  # Generate columns "target_mirna_pos" / "target_mirna_len_pos" for "Positions_Gap == 0"
  miranda_pos_1$Positions_expanded <- miranda_pos_1$Positions
  
  miranda_pos_1$target_mirna_pos <- paste(miranda_pos_1$target_mirna, paste0("pos", str_trim(miranda_pos_1$Positions)), sep="_")
  miranda_pos_1$target_mirna_len_pos <- paste(miranda_pos_1$target_mirna_len, paste0("pos", str_trim(miranda_pos_1$Positions)), sep="_")

  # Generate columns "target_mirna_pos" / "target_mirna_len_pos" for "Positions_Gap > 0"
  miranda_pos_not1_position <- ""
  for (i in 1:nrow(miranda_pos_not1)) {
    miranda_pos_not1_position <- c(miranda_pos_not1_position, unlist(str_split(miranda_pos_not1$Positions," ")[[i]]))
  }  
  miranda_pos_not1_target <- ""
  for (i in 1:nrow(miranda_pos_not1)) {
    miranda_pos_not1_target <- c(miranda_pos_not1_target, 
                                 unlist(paste0(miranda_pos_not1$target_mirna_len[i], "_pos", unlist(str_split(miranda_pos_not1$Positions," ")[[i]])))
                                 )
  }
 
  miranda_pos_not1_expanded <- data.frame(Positions_expanded = miranda_pos_not1_position[-1], 
                                          target_mirna_len_pos = miranda_pos_not1_target[-1])
  miranda_pos_not1_expanded$target_mirna_len <- sub("_pos.+", "", miranda_pos_not1_expanded$target_mirna_len_pos)

  miranda_pos_not1_expanded <- miranda_pos_not1_expanded %>% merge(miranda_pos_not1,by="target_mirna_len") # data frame with 15 columns
  miranda_pos_not1_expanded$target_mirna_pos <- paste(miranda_pos_not1_expanded$target_mirna, 
                                                      paste0("pos",miranda_pos_not1_expanded$Positions_expanded), sep="_")
}
head(miranda_pos_1)
head(miranda_pos_not1_expanded)

# miRanda: combine data frames "miranda_pos_1" and "miranda_pos_not1_expanded"
{
  miranda_positioned <- bind_rows(miranda_pos_1,miranda_pos_not1_expanded)
  miranda_positioned <- miranda_positioned %>% relocate(Positions_expanded, .after=Positions)
  miranda_positioned$Positions_expanded <- as.numeric(miranda_positioned$Positions_expanded)
} 

write.csv(miranda_positioned, paste0("Tco_Cody.3utr.fa.miranda_with_strict.hit_positioned.csv"),row.names=FALSE)
write.csv(miranda_positioned, paste0("Hho_Cody.3utr.fa.miranda_with_strict.hit_positioned.csv"),row.names=FALSE)
write.csv(miranda_positioned, paste0("Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.miranda_with_strict.hit_positioned.csv"),row.names=FALSE)
write.csv(miranda_positioned, paste0("Dme_3_prime_UTR_dedup.fasta.miranda_with_strict.hit_positioned.csv"),row.names=FALSE)
####

# Read expanded files of miRanda
miranda_positioned <- read.csv("Tco_Cody.3utr.fa.miranda_with_strict.hit_positioned.csv")
miranda_positioned <- read.csv("Hho_Cody.3utr.fa.miranda_with_strict.hit_positioned.csv")
miranda_positioned <- read.csv("Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.miranda_with_strict.hit_positioned.csv")
miranda_positioned <- read.csv("Dme_3_prime_UTR_dedup.fasta.miranda_with_strict.hit_positioned.csv")
head(miranda_positioned)


# Validate corrected length
validate_length2 <- function(RNAhybrid_217, miranda) {
  RNAhybrid_target_unique <- RNAhybrid_217_positioned[!duplicated(RNAhybrid_217_positioned$target_id),]
  miranda_target_unique <- miranda_positioned[!duplicated(miranda_positioned$target),]
  
  RNAhybrid_miranda_target_unique <- merge(RNAhybrid_target_unique,miranda_target_unique,by.x="target_id",by.y="target")
  head(RNAhybrid_miranda_target_unique)
  print(RNAhybrid_miranda_target_unique$target_length_corrected == RNAhybrid_miranda_target_unique$target_Len)
  
  RNAhybrid_miranda_target_unique <- 
    RNAhybrid_miranda_target_unique[!RNAhybrid_miranda_target_unique$target_length_corrected == RNAhybrid_miranda_target_unique$target_Len, ]
  return(RNAhybrid_miranda_target_unique) 
}
RNAhybrid_miranda_target_positioned <- validate_length2(RNAhybrid_217_positioned, miranda_positioned)

# Compare the two data frames:RNAhybrid and miRanda by "target_mirna_pos" 
# (perfect match of "Positions") (position_match_range <- 0)
{
  RNAhybrid_217_miranda_positioned_structured <- merge(RNAhybrid_217_positioned_structured, miranda_positioned, by="target_mirna_pos")
  # by "target_mirna_len_pos"
  summary(RNAhybrid_217_miranda_positioned_structured$target_mirna_len_pos.x == RNAhybrid_217_miranda_positioned_structured$target_mirna_len_pos.y)
  
  head(RNAhybrid_217_miranda_positioned_structured)
}
write.csv(RNAhybrid_217_miranda_positioned_structured,paste0(org, "_RNAhybrid_2-7_12-17_miranda_positioned_structured_ranged0.csv"), row.names=FALSE)
####

# Compare the two data frames:RNAhybrid and miRanda by "target_mirna"
RNAhybrid_217_miranda_positioned_structured_expanded <- merge(RNAhybrid_217_positioned_structured,miranda_positioned,by="target_mirna")
write.csv(RNAhybrid_217_miranda_positioned_structured_expanded, paste0(org, "_RNAhybrid_2-7_12-17_miranda_positioned_structured_expanded.csv"), row.names=FALSE)
####

RNAhybrid_217_miranda_positioned_structured_expanded <- read.csv(paste0(org,"_RNAhybrid_2-7_12-17_miranda_positioned_structured_expanded.csv"))

# Match position range
position_match_range <- 22
{
  position_corrected_min <- RNAhybrid_217_miranda_positioned_structured_expanded$position_corrected - position_match_range
  position_corrected_max <- RNAhybrid_217_miranda_positioned_structured_expanded$position_corrected + position_match_range
  # match "Positions_expanded" (miRanda) to "Positions_range" (RNAhybrid)
  RNAhybrid_217_miranda_positioned_structured_expanded_ranged <- 
    RNAhybrid_217_miranda_positioned_structured_expanded [
      (RNAhybrid_217_miranda_positioned_structured_expanded$Positions_expanded >= position_corrected_min) & 
      (RNAhybrid_217_miranda_positioned_structured_expanded$Positions_expanded <= position_corrected_max),
      ]
}
# Do not need to remove mRNA isoforms with identical 3'UTR sequences

# Do not need to remove mRNA isoforms (noisoform)

# continue on "miRNA_target_prediction_pipeline.R"

# -> step 5a. Retrieve miRNAs annotation and sequences from source file, annotate miRNAs 
# RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_5a <- annotate_mirna(RNAhybrid_217_miranda_positioned_structured_expanded_ranged)

# -> step 5b. Retrieve target mRNA gene name and sequences from source file, annotate target mRNA
# RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq <- 
#   annotate_mrna(RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_5a, position_match_range = position_match_range)

# -> step 2. Remove mRNA isoform with identical 3'UTR sequences, make miRNA and mRNA binding structure for viewing in Excel, calculate numbers of bulges/loops

# Do not need to remove mRNA isoforms with identical 3'UTR sequences
# dedup <- read.csv(paste0(org, "_dedup.csv")) # dataframe of two columns ("Gene" and "mRNA_id") containing all non-redundant mRNA transcripts
# {
#   TF_RNAhybrid_miranda_dedup <- RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq$target_id %in% dedup$mRNA_id
#   summary(TF_RNAhybrid_miranda_dedup)
#   RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq_dedup <- RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq[TF_RNAhybrid_miranda_dedup, ]
#   
#   RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq_dedup <- RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq_dedup %>% relocate(target_mirna,.before=target_id)
#   head(RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq_dedup)
# }

# -> step 4. Remove all mRNA isoforms
# Would lose information on the predicted binding position, only do this for finding out which miRNAs (Nmirna) or targets (Ntarget) are present in this dataset
# RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq_noisoform <- isoform_remover(RNAhybrid_217_miranda_positioned_structured_expanded_ranged_annotated_seq)

# -> step 6a. Sort by miRNA occurrence (optional)
# -> step 6b. Sort by mRNA occurrence (optional)

##End OF PIPELINE
