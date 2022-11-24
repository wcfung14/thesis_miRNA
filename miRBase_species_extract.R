#' @title Extract miRNA hairpin from miRBase for a species
#' @description Extract Dme miRNA hairpin from all miRBase miRNA sequences "Fasta format sequences of all miRNA hairpins" (https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz)

mirbase_extract <- function(org) {
  mirna_fasta <- seqinr::read.fasta(file = "hairpin.fa")
  message(paste0("Number of sequences found: "), sum(grepl(org, names(mirna_fasta), ignore.case=TRUE)))
  
  mirna_fasta_extract <- mirna_fasta[c(grep(org, names(mirna_fasta), ignore.case=TRUE))]
  
  seqinr::write.fasta(mirna_fasta_extract, names = names(mirna_fasta_extract), file.out = paste0(org, "_hairpin.fa"))
  message(paste0("File output as: ", org, "_hairpin.fa"))
  return()
}


# runs only when script is run by itself
if (sys.nframe() == 0){
  setwd(here::here())
  
  org = "Dme"
  setwd(paste0("../", org))
  mirbase_extract(org)
}
