#' @title Extract miRNA hairpin from miRBase for a species
#' @description Extract miRNA hairpin from all miRBase miRNA sequences "Fasta format sequences of all miRNA hairpins" (https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz)
#' For Dme

mirbase_extract <- function(org) {
  mirna_fasta <- seqinr::read.fasta(file = "hairpin.fa")
  message(paste0("Number of sequences found: "), sum(grepl(org, names(mirna_fasta), ignore.case=TRUE)))
  
  mirna_fasta_extract <- mirna_fasta[c(grep(org, names(mirna_fasta), ignore.case=TRUE))]
  
  seqinr::write.fasta(mirna_fasta_extract, names = names(mirna_fasta_extract), file.out = paste0(org, "_hairpin.fa"))
  message(paste0("File is outputed as: ", org, "_hairpin.fa",
                 "\n at path: \n", getwd()))
  return()
}


# runs only when script is run by itself
if (sys.nframe() == 0){
  org = "Dme"
  mirbase_extract(org)
}
