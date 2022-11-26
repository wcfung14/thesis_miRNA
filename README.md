# Thesis - miRNA target prediction
R script for miRNA target prediction for my thesis submitted in 2022

## script
`1-6.miRNA_target_prediction_pipeline.R` - miRNA target prediction pipeline
`1.miRNA_target_prediction_pipeline_position.R` -  miRNA target prediction by position of prediction (to correct for a bug found in dataset of RNAhybrid)
`5a.miRBase_species_extract.R` -  Extract miRNA hairpin from miRBase for a species (for Dme)
`5a.miRDeep2_annotation.R` - miRNA target prediction annotation by miRDeep2 (for Tco, Hho, Nda)
`7a.upset_Dme.R` -  Confident miRNA-target predictions (for Dme)
`7b.venn.R` - Compare miRNA-target predictions across species by Venn Diagram
`8.miRNA_prediction_summary_species_to_miRNA_formatter.R` - miRNA prediction summary formatting from each species to each miRNA (optional)
