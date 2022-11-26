#' @title Compare miRNA-target predictions across species by Venn Diagram
#' @description Compare miRNA-target predictions across species by Venn Diagram
#' Include all species

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(here::here())
setwd(paste0("../"))

# load required library
library(xlsx)
library(ggvenn)
library(UpSetR)

{
  position_match_range <- 99999
  Tco_miRNAs <- read.csv(paste("Tco/Tco_RNAhybrid_2-7_12-17_miranda_positioned_structured_expanded_ranged", position_match_range,"_annotated_seq.csv", sep=""))
  Hho_miRNAs <- read.csv(paste("Hho/Hho_RNAhybrid_2-7_12-17_miranda_positioned_structured_expanded_ranged", position_match_range,"_annotated_seq.csv", sep=""))
  Nda_miRNAs <- read.csv(paste("Nda/Nda_RNAhybrid_2-7_12-17_miranda_positioned_structured_expanded_ranged", position_match_range,"_annotated_seq.csv", sep=""))
  Dme_miRNAs <- read.xlsx(paste("Dme/Dme_miRNA_summary_UpSet_witharms.xlsx"), sheetName = "Dme_UpSet") # result from "upset_Dme.R"
  Dme_miRNAs <- read.xlsx(paste("Dme/Dme_miRNA_summary_UpSet_noarms.xlsx"), sheetName = "Dme_UpSet")# result from "upset_Dme.R"
}

dme_formatting <- function(df) {
  df <- df[df$degree == 5, "Dme_mirna_list_un"]
  df <- gsub("dme-|", "", df)
  # df <- gsub("-.p", "", df)
  gene_mirna <- gsub("miR", "mir", df)
  return(gene_mirna)
}
Dme_miRNAs_gene_mirna <- dme_formatting(Dme_miRNAs)

formatting <- function(df){
  df <- df[, c("Gene", "miRNA_id", "miRNA_annotation")]
  df <- df[!grepl(c("Novel|no hits|NULL|MIR+"), df$miRNA_annotation),]
  df$arm <- stringr::str_sub(df$miRNA_id, start=-3)
  gene_mirna <- paste(paste0(df$miRNA_annotation, df$arm),
                      df$Gene, sep="_")
  return(gene_mirna)
}
Tco_miRNAs_gene_mirna <- formatting(Tco_miRNAs)
Hho_miRNAs_gene_mirna <- formatting(Hho_miRNAs)
Nda_miRNAs_gene_mirna <- formatting(Nda_miRNAs)

mirna_list <- list(Dme = Dme_miRNAs_gene_mirna, Tco = Tco_miRNAs_gene_mirna, 
                   Hho = Hho_miRNAs_gene_mirna, Nda = Nda_miRNAs_gene_mirna)

# Draw venn diagram
draw_venn <- function(mirna_list = mirna_list) {
  grid.newpage()
  ggvenn(mirna_list)
}
draw_venn(mirna_list)

# Draw upset plot
upset_plot <- upset(fromList(mirna_list), nsets=length(mirna_list), nintersect=NA, order.by = c("degree"))
upset_plot

# Make result dataframe (with UpsetR)
# result format: table with the first column representing miRNA-target predictions 
# and subsequent columns representing a dataset (algorithm) (dataset name in first row) 
# "1" = the miRNA-target prediction was predicted by the corresponding algorithms
# laste two columns: "concat" and "degree"
{
  # concat results into one string (e.g. "X1111")
  upset_matrix <- upset_plot$New_data
  upset_matrix$concat <- do.call(paste0, upset_matrix[1:ncol(upset_matrix)])
  table(upset_matrix$concat)
  
  # formatting the result by unlisting
  mirna_list_un <- unlist(mirna_list, use.names = FALSE)
  mirna_list_un <- mirna_list_un[!duplicated(mirna_list_un)]
  
  # make dataframe
  upset_df <- data.frame(mirna_list_un, upset_matrix)
}

# filter resutls (optional, easier to do in Excel instead)
{
  upset_filter <- which(upset_df$concat %in% "1111")
  upset_df_filtered <- upset_df[upset_filter, ]
}


# Retrieve list of items in each set and intersection
{
  v.table <- gplots::venn(mirna_list)
  # print(v.table)
  # unlist(v.table)
  
  names(attributes(v.table))
  attr(x = v.table, "intersections")
  attr(x = v.table, "intersections")$`Dme:Tco`
}