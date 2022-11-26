#' @title miRNA target prediction
#' @description Combine results from RNAhybrid and miRanda, remove target isoforms with identical 3'UTR seq, and filter for internal/bulge loop, return summary of results

org = "Tco" # c("Dme", "Tco", "Hho", "Nda")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(here::here())
setwd(paste0("../", org))

# load required library
library(xlsx)
library(tidyverse)
library(Biostrings)

{
  # Read result files of RNAhybrid and miRanda (.xls)
  # A bug was found (on 20210823) for results of RNAhybrid: for every 1000 bp of the 3'UTR,
  # the predicted miRNA binding sites location will be decreased by 1 bp
  # see alternative "miRNA_target_prediction_pipeline_position.R" for position correction 
  RNAhybrid_217 <- read.delim("Tco_Cody.3utr.fa.RNAhybrid_targets_seed-region_2-7_12-17.xls", stringsAsFactors=FALSE, sep="\t")
  RNAhybrid_217 <- read.delim("Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.RNAhybrid_targets_seed-region_2-7_12-17.xls", stringsAsFactors=FALSE, sep="\t")
  RNAhybrid_217 <- read.delim("Dme_3_prime_UTR_dedup.fasta.RNAhybrid_targets_seed-region_2-7_12-17.xls", stringsAsFactors=FALSE, sep="\t")
  head(RNAhybrid_217)
  miranda <- read.delim("Tco_Cody.3utr.fa.miranda_with_strict.hit.xls", stringsAsFactors=FALSE, sep="\t")
  miranda <- read.delim("Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta.miranda_with_strict.hit.xls", stringsAsFactors=FALSE, sep="\t")
  miranda <- read.delim("Dme_3_prime_UTR_dedup.fasta.miranda_with_strict.hit.xls", stringsAsFactors=FALSE, sep="\t")
  head(miranda)
}

# 1. Concatenate "target_id" and "miRNA_id"
concatenate_target_miRNA <- function(RNAhybrid_217, miranda) {
  RNAhybrid_217$target_mirna <- paste(RNAhybrid_217$target_id,RNAhybrid_217$miRNA_id,sep="_")
  miranda$target_mirna <- paste(miranda$target,sub("^>>", "", miranda$miRNAs),sep="_")

  ##Compare the two data frames
  TF_RNAhybrid_miranda <- RNAhybrid_217$target_mirna %in% miranda$target_mirna
  message("Number of common predictions (TRUE) by RNAhybrid and miRanda: ")
  print(summary(TF_RNAhybrid_miranda))
  
  RNAhybrid_217_miranda <- RNAhybrid_217[TF_RNAhybrid_miranda,]
  return(RNAhybrid_217_miranda)
}
RNAhybrid_217_miranda <- concatenate_target_miRNA(RNAhybrid_217, miranda)
#write.csv(RNAhybrid_217_miranda,file=paste0(org,"_RNAhybrid_2-7_12-17_miranda.csv"),row.names=FALSE)
  
# 2. Remove mRNA isoform with identical 3'UTR sequences, make miRNA and mRNA binding structure for viewing in Excel, calculate numbers of bulges/loops
dedup <- read.csv(paste0(org, "_dedup.csv")) # dataframe of two columns ("Gene" and "mRNA_id") containing all non-redundant mRNA transcripts
view_structure <- function(RNAhybrid_217_miranda, dedup) {
  library(dplyr)
  library(stringr)
  
  # Remove mRNA isoforms with identical 3'UTR sequences (have to check isoforms alignment manually beforehand)
  TF_RNAhybrid_miranda_dedup <- RNAhybrid_217_miranda$target_id %in% dedup$mRNA_id
  message("Number of predictions (TRUE) by RNAhybrid and miRanda after identical 3'UTR removal: ")
  print(summary(TF_RNAhybrid_miranda_dedup))
  RNAhybrid_217_miranda_dedup <- RNAhybrid_217_miranda[TF_RNAhybrid_miranda_dedup, ]
  RNAhybrid_217_miranda_dedup <- RNAhybrid_217_miranda_dedup %>% relocate(target_mirna) # move "target_mirna" column to the front
  
  # Concatenate target and miRNA alignment for viewing binding structure
  miRNA_structure <- paste(RNAhybrid_217_miranda_dedup$target_aln, RNAhybrid_217_miranda_dedup$target_hit, RNAhybrid_217_miranda_dedup$miRNA_hit, RNAhybrid_217_miranda_dedup$miRNA_aln, sep="\n")
  RNAhybrid_217_miranda_dedup_structure <- cbind(RNAhybrid_217_miranda_dedup,miRNA_structure)
  message("To view miRNA and target binding structure:")
  message('In Excel, change font to "Consolas" for the column "Structure" and wrap text to view alignment')
  
  # Calculate the number of internal bulges AND loops in miRNA and target binding structure
  # Sorry I have not figured out how to filter only bulges OR internal loops OR perfect internal loops
  RNAhybrid_217_miranda_dedup_structure$internal_bulge_loop <- 
    str_count(str_squish(str_trim(RNAhybrid_217_miranda_dedup$miRNA_hit))," ")
  head(RNAhybrid_217_miranda_dedup_structure)
  message("Distribution of number of internal bulge and loops in miRNA and target binding structure: ")
  print(count(RNAhybrid_217_miranda_dedup_structure, internal_bulge_loop))
  
  # Initialize columns for manual filtering on miRNA and target binding structure in Excel
  RNAhybrid_217_miranda_dedup_structure$internal_loop <- "" 
  RNAhybrid_217_miranda_dedup_structure$perfect_internal_loop <- ""
  RNAhybrid_217_miranda_dedup_structure$structure_okay <- ""
  message('In Excel, filter column "internal loop"=1, check for typical mRNA target and miRNA binding structure (internal bulge/loop between 12 and 17)')
  message('Input 1 in column "structure_okay" if structure is okay for further filtering')
  message('rename the file to "...structured" after manual filtering')
  message(paste0('Result is outputed as ', org, "_RNAhybrid_2-7_12-17_miranda_dedup_structure.xlsx", 
                 "\n at path: \n", getwd()))
  
  write.xlsx(RNAhybrid_217_miranda_dedup_structure,file = paste0(org, "_RNAhybrid_2-7_12-17_miranda_dedup_structure.xlsx"), row.names=FALSE)
  return(RNAhybrid_217_miranda_dedup_structure)
}
RNAhybrid_217_miranda_dedup_structure <- view_structure(RNAhybrid_217_miranda, dedup)

# 3. Filter for typical miRNA and target binding structure after manual checking ("structure_okay"=1)
# File was renamed to "...structured" after manual filtering
RNAhybrid_217_miranda_dedup_structured <- read.xlsx(paste0(org,"_RNAhybrid_2-7_12-17_miranda_dedup_structured.xlsx"), sheetIndex=1)

loop_filter <- function(RNAhybrid_217_miranda_dedup_structured, loop_freq_plot = TRUE){

  #table(RNAhybrid_217_miranda_dedup_structured$internal_bulge_loop)
  #table(RNAhybrid_217_miranda_dedup_structured$internal_loop)
  #table(RNAhybrid_217_miranda_dedup_structured$perfect_internal_loop)
  #table(RNAhybrid_217_miranda_dedup_structured$structure_okay)

  # Plot internal_bulge_loop frequencies
  if (loop_freq_plot) {
    barplot_x <- RNAhybrid_217_miranda_dedup_structured$internal_bulge_loop
    ylim <- c(0, 500 * (ceiling(max(table(barplot_x)) / 500)))
    ylim <- c(0, 1.1 * max(table(barplot_x)))
    xx <- barplot(table(barplot_x), width=0.85, ylim=ylim, main=org, cex.axis=0.8, xlab="Number of internal/bulge loops", ylab="Frequency")
    text(x=xx, y=table(barplot_x), label=table(barplot_x), pos=3,cex=0.8)
  }
  
  # Filtering parameters
  INTERNAL_BULGE_LOOP_QUERY <- 0:6
  INTERNAL_LOOP_QUERY <- c(NA,"", "0", "1")
  PERFECT_INTERNAL_LOOP_QUERY <- c(NA, "", "0", "1")
  STRUCTURE_QUERY <- c("1")
  LOOP <- "IBLall"

  # filter by parameters
  RNAhybrid_217_miranda_dedup_structured_looped <- RNAhybrid_217_miranda_dedup_structured [
    RNAhybrid_217_miranda_dedup_structured$internal_bulge_loop %in% INTERNAL_BULGE_LOOP_QUERY 
    & RNAhybrid_217_miranda_dedup_structured$internal_loop %in% INTERNAL_LOOP_QUERY 
    & RNAhybrid_217_miranda_dedup_structured$perfect_internal_loop %in% PERFECT_INTERNAL_LOOP_QUERY 
    & RNAhybrid_217_miranda_dedup_structured$structure_okay %in% STRUCTURE_QUERY,
  ]
  message(paste0("number of miRNA-target prediction satisfying the parameters: "), 
          nrow(RNAhybrid_217_miranda_dedup_structured_looped))
  return(RNAhybrid_217_miranda_dedup_structured_looped)
}
RNAhybrid_217_miranda_dedup_structured_looped <- loop_filter(RNAhybrid_217_miranda_dedup_structured, loop_freq_plot = TRUE)

# 4. Remove all mRNA isoforms (mRNA isoforms only needed for miRNA-target prediction, not needed for subsequent molecular validation)
isoform_remover <- function(RNAhybrid_217_miranda_dedup_structured_looped) {
  if (org == "Dme"){
    # Remove all mRNA isoforms - remove isoform indexes for Dme
    RNAhybrid_217_miranda_dedup_structured_looped$target_id_de <- 
      sub(":.+","",RNAhybrid_217_miranda_dedup_structured_looped$target_id)
  } else {
    # Remove all mRNA isoforms - remove isoform indexes for other arthropods (Tco, Hho and Nda)
    RNAhybrid_217_miranda_dedup_structured_looped$target_id_de <- 
      sub("_i.+","",RNAhybrid_217_miranda_dedup_structured_looped$target_id)
  }
  # Remove all mRNA isoforms - remove duplicated isoforms
  target_noisoform_mirna <- 
    paste(RNAhybrid_217_miranda_dedup_structured_looped$target_id_de, RNAhybrid_217_miranda_dedup_structured_looped$miRNA_id, sep="_")
  
  TF_target_noisoform_mirna <- duplicated(target_noisoform_mirna)
  message('Number of predicted miRNA-target predictions ("TRUE") after duplicated predictions removal: ')
  print(summary(TF_target_noisoform_mirna))
  RNAhybrid_217_miranda_dedup_structured_looped_noisoform <- RNAhybrid_217_miranda_dedup_structured_looped[!TF_target_noisoform_mirna,]
  return(RNAhybrid_217_miranda_dedup_structured_looped_noisoform)
}
RNAhybrid_217_miranda_dedup_structured_looped_noisoform <- isoform_remover(RNAhybrid_217_miranda_dedup_structured_looped)

# 5a. Retrieve miRNAs annotation and sequences from source file, annotate miRNAs
annotate_mirna <- function(RNAhybrid_217_miranda_dedup_structured_looped_noisoform) {
  if (org == "Dme") {
    # Retrieve miRNAs annotation and sequences from source file for Dme
    miRNA_fasta_df <- readRNAStringSet(paste0(org, "_mature.fa")) # parsed with "miRBase_species_extract.R"
    miRNA_fasta_df_annotated <- data.frame(miRNA_fasta_names = names(miRNA_fasta_df), miRNA_fasta_seq = paste(miRNA_fasta_df))
    # Annotate miRNAs for Dme
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq <- 
      RNAhybrid_217_miranda_dedup_structured_looped_noisoform %>%
      merge(miRNA_fasta_df_annotated, by.x="miRNA_id", by.y="miRNA_fasta_names") %>%
      relocate(c("miRNA_fasta_seq"), .after=miRNA_id)
    
  } else {
    # Retrieve miRNAs annotation and sequences from source file for other arthropods (Tco, Hho and Nda)
    miRNA_fasta_df_annotated <- read.xlsx(paste0(org, "_miRDeep2.mature_annotated.xlsx"), sheetIndex=1) # parsed with "miRDeep2_annotation.R"
    miRNA_fasta_df_annotated <- miRNA_fasta_df_annotated[!duplicated(miRNA_fasta_df_annotated$miRNA_fasta_names), ]
    
    # Annotate miRNAs for other arthropods (Tco, Hho and Nda)
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq <-
      RNAhybrid_217_miranda_dedup_structured_looped_noisoform %>%
      merge(miRNA_fasta_df_annotated, by.x="miRNA_id", by.y="miRNA_fasta_names") %>%
      relocate(c("miRNA_fasta_seq","miRNA_annotation"), .after=miRNA_id)
  }
  return(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq)
}
RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5a <- annotate_mirna(RNAhybrid_217_miranda_dedup_structured_looped_noisoform)

# 5b. Retrieve target mRNA gene name and sequences from source file, annotate target mRNA
annotate_mrna <- function(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5a, loop = LOOP) {
  {
    gene_name <- read.csv(paste0(org, "_query.csv")) # dataframe of containing all mRNA transcripts for each gene
    gene_name <- gene_name[, c("mRNA_id","Gene")] # Just need to transcipt ID and gene name for gene name annotation
  }
  
  # Retrieve target mRNA sequences from source file
  UTR_bed <- c(
    "Tco" = "Tco_10samples.Trinity-GG.fasta.3_prime_UTR.bed.fa",
    "Hho" = "Hho_10samples.Trinity-GG.fasta.3_prime_UTR.bed.fa",
    "Nda" = "Nda_Cody.Trinity.fasta.3_prime_UTR.bed.fasta",
    "Dme" = "Dme_3_prime_UTR_dedup.fasta"
  )
  target_fasta <- readDNAStringSet(UTR_bed[[org]])
  target_fasta_df <- data.frame(target_fasta_names = names(target_fasta), target_fasta_seq = paste(target_fasta))

  if (org == "Dme") {
    # Retrieve target mRNA sequences from source file for Dme
    target_fasta_df$target_fasta_names <- sub(":.+","",target_fasta_df$target_fasta_names)
    target_fasta_df <- target_fasta_df[!duplicated(target_fasta_df$target_fasta_names),]

    # Annotate target mRNA for Dme
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5ab <-
      RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5a %>%
      merge(gene_name, by.x="target_id_de", by.y="mRNA_id",all=FALSE) %>%
      merge(target_fasta_df, by.x="target_id_de", by.y="target_fasta_names", all=FALSE) %>%
      relocate(c("Gene","target_fasta_seq"),.after=target_id) %>%
      relocate("target_mirna", .before=target_id) %>%
      relocate("target_id_de", .after=target_id)
  
  } else {
    # Retrieve target mRNA sequences from source file for other arthropods (Tco, Hho and Nda)
    target_fasta_df <- target_fasta_df[!duplicated(target_fasta_df$target_fasta_names), ]

    # Annotate target mRNA for other arthropods (Tco, Hho and Nda)
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5ab <- 
      RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5a %>%
      merge(gene_name, by.x="target_id_de", by.y="mRNA_id",all=FALSE) %>%
      merge(target_fasta_df, by.x="target_id", by.y="target_fasta_names", all=FALSE) %>%
      relocate(c("Gene","target_fasta_seq"),.after=target_id) %>%
      relocate("target_mirna", .before=target_id)
  }
  message(paste0('Result is outputed as ', org, "_RNAhybrid_2-7_12-17_miranda_dedup_structured_", loop, 
                 "_noisoform_annotated_seq.csv", "\n at path: \n", getwd()))
  write.csv(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5ab, 
            file = paste0(org, "_RNAhybrid_2-7_12-17_miranda_dedup_structured_", LOOP, "_noisoform_annotated_seq.csv"), row.names=FALSE)
  return(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_5ab)
}
RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq <- 
  annotate_mrna(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_seq_5a, loop = LOOP)


# 6a. Sort by miRNA occurrence (optional)
sort_by_miRNA <- function(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq) {
  RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Nmirna <- 
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq %>% 
    count(miRNA_id, sort=TRUE) %>% 
    merge(miRNA_fasta_df_annotated, by.x="miRNA_id", by.y="miRNA_fasta_names") %>%
    arrange(desc(n), miRNA_id)
  
  message(paste0('Result is outputed as ', org, "_RNAhybrid_2-7_12-17_miranda_dedup_structured_", loop, 
                 "_noisoform_annotated_seq_Nmirna.csv", "\n at path: \n", getwd()))
  write.csv(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Nmirna, 
            file = paste0(org,"_RNAhybrid_2-7_12-17_miranda_dedup_structured_", LOOP, "_noisoform_annotated_seq_Nmirna.csv"), row.names=FALSE)
}
sort_by_miRNA(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq)
  
# 6b. Sort by mRNA occurrence (optional)
sort_by_mRNA <- function(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq) {
 
   RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget <- 
     RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq %>% 
     count(target_id, sort=TRUE) 

  if (org == "Dme") {
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget$target_id_de <- 
      sub(":.+","", RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget$target_id)
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget <- 
      RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget %>% 
      merge(gene_name, by.x="target_id_de", by.y="mRNA_id",all=FALSE) %>% 
      merge(target_fasta_df, by.x="target_id_de", by.y="target_fasta_names") %>% 
      relocate("target_id_de", .after=target_id) %>%
      arrange(desc(n), target_id)
  } else {
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget$target_id_de <- 
      sub("_i.+","", RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget$target_id)
    RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget <- 
      RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget %>% 
      merge(gene_name, by.x="target_id_de", by.y="mRNA_id",all=FALSE) %>%
      merge(target_fasta_df, by.x="target_id", by.y="target_fasta_names") %>% 
      subset(select=-target_id_de) %>%
      arrange(desc(n), target_id)
  }
   message(paste0('Result is outputed as ', org, "_RNAhybrid_2-7_12-17_miranda_dedup_structured_", loop, 
                  "_noisoform_annotated_seq_Ntarget.csv", "\n at path: \n", getwd()))
  write.csv(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq_Ntarget, 
            file = paste0(org,"_RNAhybrid_2-7_12-17_miranda_dedup_structured_", LOOP, "_noisoform_annotated_seq_Ntarget.csv"), row.names=FALSE)
}
sort_by_mRNA(RNAhybrid_217_miranda_dedup_structured_looped_noisoform_annotated_seq)

## End OF PIPELINE

# Data analysis:
# Check target of a miRNA -> "internal_bulge_loop_filter.R" (Not used)
# Confident prediction for Dme -> "upset_Dme.R" (7a.)
# Compare miRNAs across species -> "venn.R" (7b.)
# miRNA prediction summary formatting from each species to each miRNA (optional, not very useful) (8.)



# Compare results of R_pipeline and Galaxy in comparing RNAhybrid and miRanda data sets
{
  galaxy <- read.csv("Galaxy32-Compare_two_Datasets_on_data_17_and_data_31.tabular.csv", header=FALSE, fileEncoding="UTF-8-BOM")
  RNAhybrid_217_dedup <- read.csv("hho_217_dedup.csv", )
  RNAhybrid_217_dedup <- RNAhybrid_217_dedup[,2:12]
  head(galaxy)
  head(RNAhybrid_217_dedup)
  
  identical(galaxy, RNAhybrid_217_dedup)
  TF_galRNAhybrid <- galaxy[, 12] %in% RNAhybrid_217_dedup[, 12]
  summary(TF_galRNAhybrid)
  summary(arsenal::comparedf(galaxy,RNAhybrid_217_dedup))
}
