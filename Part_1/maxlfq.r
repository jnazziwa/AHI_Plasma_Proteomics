library(data.table)
library(tidyverse)
library(iq)

setwd("~/PhD/Projects/proteomics/Preprocessing_IQ")

raw <- data.table::fread("20220201_073608_DEPplasmaFAIMS-HIV-Durban+IAVI_CombinedLibrarySearch_210608_Report_iq.xls")
saveRDS(raw, file = "raw_iqDEP.rds")
#raw<- readRDS("raw_neat.rds") # use this incase i need to re-run the same mssats report

selected <- raw$F.ExcludedFromQuantification == "FALSE" & 
  raw$F.FrgLossType == "noloss" &
  (is.na(raw$PG.Qvalue) | raw$PG.Qvalue <= 0.01) &
  (is.na(raw$EG.Qvalue) | raw$EG.Qvalue <= 0.01)

raw <- raw[selected,]

hist(log2(raw[, "F.PeakArea"]), 100, main = "Histogram of log2 intensities", 
     col = "steelblue", border = "steelblue", freq = FALSE)

## PREPROCESS
sample_id  <- "R.FileName" 

secondary_id <- c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType")

norm_data <- iq::preprocess(raw, 
                            sample_id  = sample_id, 
                            secondary_id = secondary_id,
                            median_normalization = FALSE)

#Create protein list
protein_list <- iq::create_protein_list(norm_data)
iq::plot_protein(protein_list$"P12821;P12821-2", main = "Protein P12821;P12821-2", split = NULL) 
iq::plot_protein(protein_list$"P12821;P12821-2", main = "Protein P12821;P12821-2", cex = 0.5) 

x = names(protein_list)

par(mfrow = c(2, 2))
for (i in protein_list) {
  iq::plot_protein(i, main = paste("protein",i), cex = 0.4)
}



#create protein table
result <- iq::create_protein_table(protein_list)


iq::plot_protein(rbind(protein_list$"P12821;P12821-2", 
                       MaxLFQ = iq::maxLFQ(protein_list$"P12821;P12821-2")$estimate), 
                 main = "MaxLFQ quantification of  P12821;P12821-2", 
                 col = c(rep("gray", nrow(protein_list$"P12821;P12821-2")), "green"), 
                 split = NULL)  



#extract data
extra_names <- iq::extract_annotation(rownames(result$estimate), 
                                      raw, 
                                      annotation_columns = c("PG.Genes", "PG.ProteinNames"))

write.table(cbind(Protein = rownames(result$estimate),
                  extra_names[, c("PG.Genes", "PG.ProteinNames")],
                  result$estimate, 
                  annotation = result$annotation),
            "output-maxLFQ-annotation.txt", sep = "\t", row.names = FALSE)

## change lab IDs to sample ids
DEP_annotation <- readr::read_tsv("annotation_depletedplasma.txt")   
DEP_annotation1 <- tidyr::unite(DEP_annotation, "sampleID", BioReplicate:Condition) %>% dplyr::select(-...1)


DEP_long_df <- readr::read_tsv("output-maxLFQ-annotation.txt")
DEP_long_df <-  dplyr::select(DEP_long_df, -annotation)
DEP_wide_df <- DEP_long_df %>% 
  tidyr::pivot_longer(cols = 4:159, names_to = "Run", values_to = "intensity") %>% 
dplyr::left_join(DEP_annotation1, by = "Run") %>% 
  dplyr::select(-Run)

#tidyr::pivot_wider(names_from = sampleID, values_from = intensity)


#####-------------------------------------------------------------------------------------------------------------
##neat


rawN <- fread("20220201_110528_ACUTEHIV_Report_iq-xls.xls")
saveRDS(rawN, file = "raw_iqNEAT.rds")
#raw<- readRDS("raw_iqNEAT.rds") # use this incase i need to re-run the same mssats report

selectedN <- rawN$F.ExcludedFromQuantification == "FALSE" & 
  rawN$F.FrgLossType == "noloss" &
  (is.na(rawN$PG.Qvalue) | rawN$PG.Qvalue <= 0.01) &
  (is.na(rawN$EG.Qvalue) | rawN$EG.Qvalue <= 0.01)

rawN <- rawN[selectedN,]


## PREPROCESS
sample_id  <- "R.FileName" 

secondary_id <- c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType")

norm_dataN <- iq::preprocess(rawN, 
                            sample_id  = sample_id, 
                            secondary_id = secondary_id,
                            median_normalization = FALSE)

#Create protein list
protein_listN <- iq::create_protein_list(norm_dataN)

iq::plot_protein(protein_listN$P05543, main = "Protein P05543", split = NULL) 
iq::plot_protein(protein_listN$P05543, main = "Protein P05543", cex = 0.8) 


#create protein table
resultN <- iq::create_protein_table(protein_listN)


iq::plot_protein(rbind(protein_listN$P05543, 
                       MaxLFQ = iq::maxLFQ(protein_listN$P05543)$estimate), 
                 main = "MaxLFQ quantification of  P05543-NEAT", 
                 col = c(rep("gray", nrow(protein_listN$P05543)), "green"), 
                 split = NULL)  



#extract data
extra_namesN <- iq::extract_annotation(rownames(resultN$estimate), 
                                      rawN, 
                                      annotation_columns = c("PG.Genes", "PG.ProteinNames"))

write.table(cbind(Protein = rownames(resultN$estimate),
                  extra_namesN[, c("PG.Genes", "PG.ProteinNames")],
                  resultN$estimate, 
                  annotation = resultN$annotation),
            "output-maxLFQ_NEAT-annotation.txt", sep = "\t", row.names = FALSE)


## change lab IDs to sample ids
NEAT_annotation <- readr::read_tsv("annotation_NonDep_2.txt")   
NEAT_annotation1 <- tidyr::unite(NEAT_annotation, "sampleID", BioReplicate:Condition) %>% dplyr::select(-...1)


NEAT_long_df <- readr::read_tsv("output-maxLFQ_NEAT-annotation.txt")
NEAT_long_df <-  dplyr::select(NEAT_long_df, -annotation)
NEAT_wide_df <- NEAT_long_df %>% 
  tidyr::pivot_longer(cols = 4:156, names_to = "Run", values_to = "intensity") %>% 
  dplyr::left_join(NEAT_annotation1, by = "Run") %>% 
  dplyr::select(-Run) 
  #tidyr::pivot_wider(names_from = sampleID, values_from = intensity)

fwrite(NEAT_wide_df,"neat_20220201.txt")
fwrite(DEP_wide_df,"depleted_20220201.txt")