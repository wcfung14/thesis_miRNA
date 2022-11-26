#' @title miRNA target prediction annotation by miRDeep2
#' @description Annotate miRNAs using results of miRDeep2 (blast) and MirGeneDB "Unique structural features of miRNAs" (if any)
#' For non-model arthropods only (Tco, Hho, Nda)

org = "Tco" # c("Tco", "Hho", "Nda")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(here::here())
setwd(paste0("../", org))

library(xlsx)
library(tidyverse)
library(Biostrings)

mirdeep2_annotate <- function(org){
  # Load miRDeep2 results
  miRNA_fasta <- readDNAStringSet(paste(org, "_miRDeep2.mature.fa", sep=""))
  miRNA_fasta_names <- names(miRNA_fasta)
  miRNA_fasta_names_dearm <- sub("-.p","",miRNA_fasta_names)
  miRNA_fasta_df <- data.frame(miRNA_fasta_names,miRNA_fasta_names_dearm,miRNA_fasta_seq = paste(miRNA_fasta))
  
  if (org %in% c("Tco", "Hho")) {
    #Filter for Tco and Hho (blasted in miRBase + manually checked for MirGeneDB by Kevin and Henry)
    mirna_annotation <- read.xlsx("../pbio.3000636.s032.xlsx",
                                  sheetName = paste(org," miRNAs",sep=""))
    # pbio.3000636.s032.xlsx was retrieve from
    # S4 Data. The microRNA contents, arm usage, and predicted gene targets (pbio.3000636.s032.xlsx)
    # Qu, Z., Nong, W., So, W. L., Barton-Owen, T., Li, Y., Leung, T. C., ... & Hui, J. H. (2020). Millipede genomes reveal unique adaptations during myriapod evolution. PLoS biology, 18(9), e3000636.
    
    mirna_annotation <- mirna_annotation[!mirna_annotation[,1] == paste0(org," novel miRNA"),c(1,2)] # delete header row for "novel miRNA", extract only the first two columns
    message('numbers of miRDeep2 predicted miRNAs that passed MirGeneDB criteria ("TRUE")')
    print(summary(miRNA_fasta_df$miRNA_fasta_names_dearm %in% mirna_annotation[,1]))
    
    miRNA_fasta_df_annotated <- merge(miRNA_fasta_df,mirna_annotation,by.x="miRNA_fasta_names_dearm",by.y=1, all=TRUE)
  } else if (org == "Nda") {
    # Filter for Nda (blasted in miRBase only)
    mirna_annotation <- read.delim(paste(org, "_miRDeep2.hairpin.fa.annatation.xls",sep=""),
                                   stringsAsFactors=FALSE, sep="\t")
    mirna_annotation <- mirna_annotation[, c("candidate_id", "miRBase_blastn_hits", "X.miRNA")]
    
    # number of miRNAs with good annotation
    message('numbers of miRNAs by blast in miRBase')
    sum(!grepl(c("NULL|MIR+"), mirna_annotation$candidate_id)) # need to manually check the list too
    
    message('numbers of miRDeep2 predicted miRNAs also annotated by blast in miRBase ("TRUE")')
    summary(miRNA_fasta_df$miRNA_fasta_names_dearm %in% mirna_annotation$X.miRNA)
    
    miRNA_fasta_df_annotated <- merge(miRNA_fasta_df,mirna_annotation,by.x="miRNA_fasta_names_dearm",by.y="X.miRNA", all=TRUE)
  }  
  
  miRNA_fasta_df_annotated[is.na(miRNA_fasta_df_annotated)] <- "no hits"
  miRNA_fasta_df_annotated <- subset(miRNA_fasta_df_annotated, select=-miRNA_fasta_names_dearm)

  write.xlsx(miRNA_fasta_df_annotated,paste0(org,"_miRDeep2.mature_annotated.xlsx"), row.names=FALSE)
  message(paste0("File is outputed as: ", org, "_miRDeep2.mature_annotated.xlsx",
                 "\n at path: \n", getwd()))
  return()
}


# runs only when script is run by itself
if (sys.nframe() == 0) {
  for (org in c("Tco", "Hho", "Nda")) {
    mirdeep2_annotate(org)
  }
}
