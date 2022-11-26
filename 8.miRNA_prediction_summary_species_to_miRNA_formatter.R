#' @title miRNA prediction summary formatting from each species to each miRNA
#' @description Extract miRNA predicted targets from summary for each species ("Millipede_miRNA_date.xlsx") and output summary for each miRNA ("miRNA_IBLall_summary_Ec")
#' Reuslt format: Excel file, in each sheet: 
#' table with the first column representing genes used for prediction
#' and subsequent columns represent miRNA homologs in different species
#' "1" = the miRNA-target prediction was predicted for the species


org = "Tco" # c("Dme", "Tco", "Hho", "Nda")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(here::here())
setwd(paste0("../", org))
miRNA.xlsx <- "Millipede_miRNA_20221116.xlsx"

# Load data files and select annotated miRNAs results
library(xlsx)
load_data <- function (org) {
  summary <- read.xlsx(miRNA.xlsx, sheetIndex = paste0("Summary_", org))
 
  # retaining the first two columns, remove columns with no annotation ("NA")
  summary_annotated <- summary[,c(1:2, which(!is.na(summary[1,])))]
  return(summary_annotated)
}
{
  Tco_summary_annotated <- load_data("Tco")
  Hho_summary_annotated <- load_data("Hho")
  Nda_summary_annotated <- load_data("Nda")
  Dme_summary_annotated <- load_data("Dme")
}
{
  colnames(Tco_summary_annotated) <- sub("Scaffold","Tco_Scaffold",colnames(Tco_summary_annotated))
  colnames(Hho_summary_annotated) <- sub("SczTNLB","Hho_SczTNLB",colnames(Hho_summary_annotated))
  colnames(Nda_summary_annotated) <- sub("Scf7vLx","Nda_Scf7vLx",colnames(Nda_summary_annotated))
  colnames(Dme_summary_annotated) <- sub("dme","Dme",colnames(Dme_summary_annotated))
}
head(Tco_summary_annotated)

# Filter for internal/bulge loops, select the range for the data sets of the specified loops
loop <- "IBLall" 
IBL_filter <- function(summary_annotated) {
  loop_index <- c(grep("^n",summary_annotated[,1]),nrow(summary_annotated))
  loop_filter <- grep(loop,summary_annotated[,1])
  summary_annotated_looped <- summary_annotated[c(1:2,loop_filter:min(loop_index[loop_filter<loop_index])),]
  return(summary_annotated_looped)
}
{
  Tco_summary_annotated_looped <- IBL_filter(Tco_summary_annotated)
  Hho_summary_annotated_looped <- IBL_filter(Hho_summary_annotated)
  Nda_summary_annotated_looped <- IBL_filter(Nda_summary_annotated)
  Dme_summary_annotated_looped <- IBL_filter(Dme_summary_annotated)
}
head(Tco_summary_annotated_looped)

# miRNA query
mirna_query <- "mir-281"

# Ecdysteroids-related (in previous papers) miRNAs # "Ec"
mirna_list <- c("bantam","let-7","mir-315","mir-2", "mir-263", "mir-9", "mir-190", "mir-252", "mir-1", "mir-277", "mir-34", "mir-210", "mir-7", "mir-965", "mir-219", "mir-281", "mir-12", "mir-33")
                # "mir-8", "mir-275", "mir-281", "mir-282", "mir-12", "mir-33", "mir-10", "mir-125", "mir-317"
# Non-ecdysteroids-related miRNAs # "nonEc"
# mirna_list <- c("mir-315", "mir-153", "mir-190", "mir-2788", "mir-307", "mir-2001", "mir-1", "mir-210", "mir-124", "mir-7", "mir-71", "mir-965", "mir-993", "mir-219") 

library(dplyr)
for (mirna_query in mirna_list) {
  
  Tco_summary_annotated_looped_query <- if(length(which(Tco_summary_annotated_looped[2,] %in% mirna_query))){
    Tco_summary_annotated_looped[, c(3,which(Tco_summary_annotated_looped[2,] %in% mirna_query))]
  } else {NA}  
  Hho_summary_annotated_looped_query <- if(length(which(Hho_summary_annotated_looped[2,] %in% mirna_query))){
    Hho_summary_annotated_looped[, c(3,which(Hho_summary_annotated_looped[2,] %in% mirna_query))]
  } else {NA}
  Nda_summary_annotated_looped_query <- if(length(which(Nda_summary_annotated_looped[2,] %in% mirna_query))){
    Nda_summary_annotated_looped[, c(3,which(Nda_summary_annotated_looped[2,] %in% mirna_query))]
  } else {NA}
  Dme_summary_annotated_looped_query <- if(length(which(Dme_summary_annotated_looped[2,] %in% mirna_query))){
    Dme_summary_annotated_looped[, c(3,which(Dme_summary_annotated_looped[2,] %in% mirna_query))]
  } else {NA}
  
  # unsure usage?
  if(class(Tco_summary_annotated_looped_query) == "data.frame") {
    Tco_summary_annotated_looped_query$id <- 1:nrow(Tco_summary_annotated_looped_query)
  } else {
    Nda_summary_annotated_looped_query$id <- 1:nrow(Nda_summary_annotated_looped_query)
  }

  # "join" function is better than "merge" as "join" does not change the order of the data sets
  mirna_summary <- Reduce(function(x,y) merge(x,y,all=TRUE), list(Tco_summary_annotated_looped_query,Hho_summary_annotated_looped_query,Nda_summary_annotated_looped_query,Dme_summary_annotated_looped_query))
  # mirna_summary <- Reduce(function(x,y) full_join(x,y,all=TRUE,copy=TRUE), list(Tco_summary_annotated_looped_query,Hho_summary_annotated_looped_query,Nda_summary_annotated_looped_query))
  mirna_summary <- mirna_summary[order(mirna_summary$id),]
  mirna_summary <- if(sum(colnames(mirna_summary) %in% "x")) {subset(mirna_summary,select=-c(x))} else {mirna_summary}
  mirna_summary <- if(sum(colnames(mirna_summary) %in% "y")) {subset(mirna_summary,select=-c(y))} else {mirna_summary}
  mirna_summary <- if(sum(colnames(mirna_summary) %in% "id")) {subset(mirna_summary,select=-c(id))} else {mirna_summary}
  mirna_summary[is.na(mirna_summary)] <- "" 
  # Incomplete/missing transcripts as NA
  mirna_summary[grep("Phantom",mirna_summary[,1]),grep("Tco",colnames(mirna_summary))] <- NA
  mirna_summary[grep("Shadow|Shade",mirna_summary[,1]),grep("Nda",colnames(mirna_summary))] <- NA
  
  mirna_summary[-2,]
  
  write.xlsx(mirna_summary[-2,],paste0("Data/miRNA_",loop,"_summary.xlsx"),row.names=FALSE,sheetName = mirna_query,append=TRUE)  
}

# For genes that have more than one copy
# Hho E78 and HR4
mirna_summary[grep("HR4$",mirna_summary[,1]),grep("Hho",colnames(mirna_summary))] <- if(sum(as.numeric(unlist(mirna_summary[grep("HR4-1|HR4-2",mirna_summary[,1]),grep("Hho",colnames(mirna_summary))])),na.rm=TRUE)) {sum(as.numeric(unlist(mirna_summary[grep("HR4-1|HR4-2",mirna_summary[,1]),grep("Hho",colnames(mirna_summary))])),na.rm=TRUE)} else {""}
mirna_summary[grep("E78$",mirna_summary[,1]),grep("Hho",colnames(mirna_summary))] <- if(sum(as.numeric(unlist(mirna_summary[grep("E78-1|E78-2",mirna_summary[,1]),grep("Hho",colnames(mirna_summary))])),na.rm=TRUE)) {sum(as.numeric(unlist(mirna_summary[grep("E78-1|E78-2",mirna_summary[,1]),grep("Hho",colnames(mirna_summary))])),na.rm=TRUE)} else {""}
mirna_summary <- if(sum(grepl("HR4-+|E78-+",mirna_summary[,1]))) {mirna_summary[-grep("HR4-+|E78-+",mirna_summary[,1]),]} else {mirna_summary}
# Dme HR3 and beta Ftz-f1
mirna_summary[grep("HR3$",mirna_summary[,1]),grep("Dme",colnames(mirna_summary))] <- if(sum(as.numeric(unlist(mirna_summary[grep("HR3-A|HR3-F",mirna_summary[,1]),grep("Dme",colnames(mirna_summary))])),na.rm=TRUE)) {sum(as.numeric(unlist(mirna_summary[grep("HR3-A|HR3-F",mirna_summary[,1]),grep("Dme",colnames(mirna_summary))])),na.rm=TRUE)} else {""}
mirna_summary[grep("Ftz-f1$",mirna_summary[,1]),grep("Dme",colnames(mirna_summary))] <- if(sum(as.numeric(unlist(mirna_summary[grep("Ftz-f1-B|Ftz-f1-C",mirna_summary[,1]),grep("Dme",colnames(mirna_summary))])),na.rm=TRUE)) {sum(as.numeric(unlist(mirna_summary[grep("Ftz-f1-B|Ftz-f1-C",mirna_summary[,1]),grep("Dme",colnames(mirna_summary))])),na.rm=TRUE)} else {""}
mirna_summary <- if(sum(grepl("HR3-+|Ftz-f1-+",mirna_summary[,1]))) {mirna_summary[-grep("HR3-+|Ftz-f1-+",mirna_summary[,1]),]} else {mirna_summary}
