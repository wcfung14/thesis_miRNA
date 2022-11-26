#' @title Confident miRNA-target predictions for Dme
#' @description Combine results from PITA, EMBL, TargetScanFly, RNAhybrid and miRanda of Drosophila melanogaster
#' For a confidetn miRNA-target prediction

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(here::here())
setwd(paste0("../"))

# load required library
library(xlsx)

# PITA dataset formatting
{
  PITA_sites_dm3_0_0_ALL <- read.delim("PITA/PITA_sites_dm3_0_0_ALL.tab")
  PITA_sites_dm3_0_0_TOP <- read.delim("PITA/PITA_sites_dm3_0_0_TOP.tab")
  PITA_sites_dm3_3_15_ALL <- read.delim("PITA/PITA_sites_dm3_3_15_ALL.tab")
  PITA_sites_dm3_3_15_TOP <- read.delim("PITA/PITA_sites_dm3_3_15_TOP.tab")
  library(arsenal)
  summary(comparedf(PITA_sites_dm3_0_0_ALL,PITA_sites_dm3_3_15_ALL))
  head(PITA_sites_dm3_0_0_ALL)
  head(PITA_sites_dm3_0_0_TOP)
  table(PITA_sites_dm3_0_0_ALL$Seed)
  table(PITA_sites_dm3_0_0_TOP$Seed)
  
  PITA_Ec_Name <- data.frame(Name = c("nvd","spo","phm","dib","sad","shd","Cyp18a1","EcR","usp","Eip74EF","Eip75B","Eip78C","Hr46","Hr4","Hr38","Hr39"),
                             Gene = c("Nvd","Spook","Phantom","Disembodied","Shadow","Shade","CYP18A1","EcR","RXR/USP","E74","E75","E78","HR3","HR4","HR38","?Ftz-f1")
                            )
  
  PITA_sites_dm3_3_15_ALL_Ec <- merge(PITA_sites_dm3_3_15_ALL,PITA_Ec_Name,by="Name")
  unique(PITA_sites_dm3_3_15_ALL_Ec$Name.1)
} 
write.csv(PITA_sites_dm3_3_15_ALL_Ec,"PITA/PITA_sites_dm3_3_15_ALL_Ec.csv",row.names=FALSE)

# EMBL dataset formatting
{
  EMBL_StarkBrennecke2005_All <- read.xlsx("EMBL/EMBL_StarkBrennecke2005_All.xlsx",sheetIndex = 1)
  head(EMBL_StarkBrennecke2005_All)
  EMBL_Ec_Name <- data.frame(Name = c("nvd","spo","phm","dib","sad","shd","Cyp18a1","EcR","usp","Eip74EF","Eip75B","Eip78C","Hr46","Hr4","Hr38","Hr39"),
                             Gene = c("Nvd","Spook","Phantom","Disembodied","Shadow","Shade","CYP18A1","EcR","RXR/USP","E74","E75","E78","HR3","HR4","HR38","?Ftz-f1")
                            )
  
  EMBL_StarkBrennecke2005_All_Ec <- merge(EMBL_StarkBrennecke2005_All,EMBL_Ec_Name,by.x="gene.name",by.y="Name")
  EMBL_StarkBrennecke2005_All_Ec
}
write.csv(EMBL_StarkBrennecke2005_All_Ec,"EMBL/EMBL_StarkBrennecke2005_All_Ec.csv",row.names=FALSE)

# TargetScanFly dataset formatting
{
  setwd(paste0("C:/Users/",Sys.info()[["user"]],"/Desktop/Work/Millipede/miRNA")) ##find Ec genes by FlyBase ID
  
  TSF_Ec_id <- read.xlsx("Millipede_miRNA_20210910.xlsx",sheetIndex = "Dme_query")
  TSF_Ec_id <- TSF_Ec_id[!duplicated(TSF_Ec_id$Gene),1:3]
  head(TSF_Ec_id)
}
setwd(paste0("C:/Users/",Sys.info()[["user"]],"/Downloads/"))
{
  TargetScanFly_default <- read.delim("TargetScanFly/TargetScan_7.2_Predicted_Targets_Info.default_predictions.txt")
  head(TargetScanFly_default)
  ##TSF_Ec_Name <- data.frame(Name = c("nvd","spo","phm","dib","sad","shd","Cyp18a1","EcR","usp","Eip74EF","Eip75B","Eip78C","Hr46","Hr4","Hr38","Hr39"),
  ##                           Gene = c("Nvd","Spook","Phantom","Disembodied","Shadow","Shade","CYP18A1","EcR","RXR/USP","E74","E75","E78","HR3","HR4","HR38","?Ftz-f1")
  ##                         )
  TargetScanFly_default_Ec <- merge(TargetScanFly_default,TSF_Ec_id,by.x="Gene.ID",by.y="FlyBase.ID")
  TargetScanFly_default_Ec
}
write.csv(TargetScanFly_default_Ec,"TargetScanFly/TargetScan_7.2_Predicted_Targets_Info.default_predictions_Ec.csv",row.names=FALSE)
{
  setwd(paste0("C:/Users/",Sys.info()[["user"]],"/Downloads/"))
  TargetScanFly_default_context <- read.delim("Predicted_Targets_Context_Scores.default_predictions.txt")
  head(TargetScanFly_default_context)
  TSF_Ec_Name <- data.frame(Name = c("nvd","spo","phm","dib","sad","shd","Cyp18a1","EcR","usp","Eip74EF","Eip75B","Eip78C","Hr3","Hr4","Hr38","Hr39"),
                            Gene = c("Nvd","Spook","Phantom","Disembodied","Shadow","Shade","CYP18A1","EcR","RXR/USP","E74","E75","E78","HR3","HR4","HR38","?Ftz-f1")
                            )
  TargetScanFly_default_context_Ec <- merge(TargetScanFly_default_context,TSF_Ec_Name,by.x="Gene.Symbol",by.y="Name")
  TargetScanFly_default_context_Ec
}
write.csv(TargetScanFly_default_context_Ec,"TargetScanFly/TargetScan_7.2_Predicted_Targets_Context_Scores.default_predictions_Ec.csv",row.names=FALSE)
{
  TargetScanFly_conserved <- read.delim("Conserved_Family_Info.txt")
  head(TargetScanFly_conserved)
  TSF_Ec_Name <- data.frame(Name = c("nvd","spo","phm","dib","sad","shd","Cyp18a1","EcR","usp","Eip74EF","Eip75B","Eip78C","Hr46","Hr4","Hr38","Hr39"),
                            Gene = c("Nvd","Spook","Phantom","Disembodied","Shadow","Shade","CYP18A1","EcR","RXR/USP","E74","E75","E78","HR3","HR4","HR38","?Ftz-f1")
                           )
  TargetScanFly_conserved_Ec <- merge(TargetScanFly_conserved,TSF_Ec_Name,by.x="Gene.Symbol",by.y="Name")
  TargetScanFly_conserved_Ec
}
write.csv(TargetScanFly_conserved_Ec,"TargetScan_7.2_Conserved_Family_Info_Ec.csv",row.names=FALSE)
####

# Read source file
{
  # an Excel table with each column representing a dataset (algorithm) (dataset name in first row) 
  # and miRNA-target predictions in subseqeut rows
  Dme_mirna <- read.xlsx("Dme_miRNA_summary_UpSet_noarms.xlsx", sheetName = "Dme_noarms_list")
  Dme_mirna <- read.xlsx("Dme_miRNA_summary_UpSet_witharms.xlsx", sheetName = "Dme_witharms_list")
  head(Dme_mirna)
  
  Dme_mirna_list <- as.list(Dme_mirna)
  Dme_mirna_list <- lapply(Dme_mirna_list, function(x) x[!is.na(x)]) ##remove NA in each set

  attributes(Dme_mirna_list)
}

##Draw venn diagram
{
  library(ggvenn)
  tiff("Dme_miRNA_summary_ggvenn.tiff", width=7000, height=5000, units='px', res=700, compression="lzw")
  
  ggvenn(Dme_mirna_list, stroke_size=0.5, set_name_size=3)
  dev.off()
}

#UpSetR
install.packages("UpSetR")
devtools::install_github("hms-dbmi/UpSetR")
{
  library(UpSetR)
  Dme_upset <- upset(fromList(Dme_mirna_list), nsets=length(Dme_mirna_list), nintersect=NA, order.by = c("degree"))
  tiff("Dme_miRNA_summary_UpSetR.tiff", width=12000, height=5000, units='px', res=800, compression="lzw")
  Dme_upset
  dev.off()
  
  attributes(Dme_upset)
  head(Dme_upset$New_data)
  
  # Make result dataframe (with UpsetR)
  # result format: table with the first column representing miRNA-target predictions 
  # and subsequent columns representing a dataset (algorithm) (dataset name in first row) 
  # "1" = the miRNA-target prediction was predicted by the corresponding algorithms
  # laste two columns: "concat" and "degree"
  {
    # concat results into one string (e.g. "X1111")
    Dme_upset_matrix <- Dme_upset$New_data
    Dme_upset_matrix$concat <- do.call(paste0, Dme_upset_matrix[1:ncol(Dme_upset_matrix)])
    head(Dme_upset_matrix)
    table(Dme_upset_matrix$concat)
  
    # formatting the result by unlisting
    Dme_mirna_list_un <- unlist(Dme_mirna_list, use.names=FALSE)
    Dme_mirna_list_un <- Dme_mirna_list_un[!duplicated(Dme_mirna_list_un)]
    
    # make dataframe
    Dme_upset_df <- data.frame(Dme_mirna_list_un, Dme_upset_matrix)
  }
  # filter resutls (optional, easier to do in Excel instead)
  {
    upset_query <- data.frame(RNAhybrid.v2.1.2.4 = 1,
                              miRanda.v3.3a = 1,
                              TargetScanFly.v7.2 = 1,
                              TargetScanFly.v7.2.default.context = 1,
                              PITA.v6_ALL = 1,
                              PITA.v6_TOP = 1,
                              EMBL.2005 = 1,
                              DIANA.microT.CDS.v5 = 1,
                              PicTar.2006 = 1) # order of data sets are different in "with arm" and "no arm" data
    upset_query$concat <- do.call(paste0, upset_query[1:ncol(upset_query)])
    upset_query
    upset_filter <- which(Dme_upset_matrix$concat %in% upset_query$concat)
    Dme_upset_df_filtered <- Dme_upset_df[upset_filter, ]
  }
}

# write.csv(Dme_upset_df, file="Dme_UpSet.csv", row.names=FALSE)
write.table(Dme_upset_df, file="Dme_UpSet.txt", row.names=FALSE) ##prevent Excel from converting concat string into numeric; In excel: transform last column to text

# alternative result formatting: 
# table with each column representing unique algorithms combination ("concat", e.g. "10000", "11111") (first row)
# and miRNA-target predictions in subseqeut rows
# "11111" means the miRNA-target predictions were predicted by all algorithms
# order of dataset is important in interpreting the result!!!!
{
  Dme_upset_list <- list()
  for (i in unique(Dme_upset_matrix$concat)) {
    TF <- Dme_upset_matrix$concat %in% i
    Dme_mirna_list_un[TF]
    Dme_upset_list[i] <- list(Dme_mirna_list_un[TF])
  }
  # write.xlsx(Dme_mirna_list_un[TF],"Dme_UpSet.xlsx",row.names=FALSE,sheetName = i,append=TRUE)  
  
  Dme_upset_list_maxlen <- lapply(Dme_upset_list, function(lst) c(lst, rep("", max(lengths(Dme_upset_list)) - length(lst))))
  as.data.frame(Dme_upset_list_maxlen)
  write.csv(Dme_upset_list_maxlen, file="Dme_UpSet_list.csv", row.names=FALSE)
}

# Retrieve list of items in each set and intersection
{
  library(gplots)
  v.table <- venn(Dme_mirna_list)
  print(v.table)
  unlist(v.table)
  
  attr(x = v.table, "intersections")
  names(attr(x = v.table, "intersections"))
  unique(attr(x = v.table, "intersections")$"RNAhybrid.2.1.2.4:miRanda.v3.3a:TargetScanFly.7.2:PITA_sites_dm3_3_15_ALL:PITA_sites_dm3_3_15_TOP")
}
