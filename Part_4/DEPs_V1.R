# Load required libraries
library(UniprotR)
library(dplyr)
library(data.table)
library(readr)
library(biomaRt)
library(tidyr)
library(ggplot2)
library(forcats)
library(org.Hs.eg.db)
library(clusterProfiler)
library(magrittr)

################################################################################
# Human Database
################################################################################

## Set up the biomaRt connection to Ensembl (human dataset) and check atrributes to use
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
biomaRt::listAttributes(mart = ensembl)

# ------------------------------------------------------------------------------


ARS_v1 <- xlsx::read.xlsx("lm_ARSav0.xlsx", 1)
ARS_v1_sign <- ARS_v1 |> tidyr::separate(PG.Genes, c("Gene.symbol",NA), sep = ";")
ARS_v1_sign1 <- ARS_v1_sign |> dplyr::select(Gene.symbol, estimate, p)
knitr::kable(head(ARS_v1_sign1))

output_df <- pathfindR::run_pathfindR(ARS_v1_sign1)


output_dfBP <- pathfindR::run_pathfindR(ARS_v1_sign1, gene_sets = "GO-BP")
output_dfKG <- pathfindR::run_pathfindR(ARS_v1_sign1, gene_sets = "KEGG")
output_dfRE <- pathfindR::run_pathfindR(ARS_v1_sign1, gene_sets = "Reactome")
output_dfGO <- pathfindR::run_pathfindR(ARS_v1_sign1, gene_sets = "GO-All")


ARS_new <- xlsx::read.xlsx("lmer_ARSvisit.xlsx", 1)
ARS_newv10 <- ARS_new |> dplyr::filter(visit_pairwise == "v0 - v1")
ARS_newv10_sign <- ARS_newv10 |> dplyr::filter(p.global < 0.005  & p.value < 0.05)

ARS_newv20 <- ARS_new |> dplyr::filter(visit_pairwise == "v0 - v2")
ARS_newv20_sign <- ARS_newv20 |> dplyr::filter(p.global < 0.005  & p.value < 0.05)

ARS_new_sign <- ARS_new |> dplyr::filter(p.global < 0.005  & p.value < 0.05)



################################################################################
#### DIFFERENTIAL ANALYSIS
################################################################################

#Input data from LMER analysis
# Read the data
linear_df <- xlsx::read.xlsx("lmer_visit_both_adjPC12.xlsx",1)

#Remove double accessions
linear_df$uniprot <- gsub("\\..*","",linear_df$Protein)

#A. Visit 1 - 0 Analysis

# STEP 1: Select columns of interest for V1-V0
V10 <- linear_df %>%
  dplyr::select(uniprot, PG.Genes, q.global, p.global, exp,
                estimate_v1...v0_IAVI, estimate_v1...v0_Durban,
                SE_v1...v0_IAVI, SE_v1...v0_Durban,
                df_v1...v0_IAVI, df_v1...v0_Durban)



# STEP 2: Filter for significant p-values and estimates in the same direction for both cohorts
V10_sign <- V10 %>%
  dplyr::filter(p.global < 0.005 & q.global < 0.005) %>%
  dplyr::filter(purrr::map2_lgl(estimate_v1...v0_Durban, estimate_v1...v0_IAVI, ~ (.x > 0 & .y > 0) | (.x < 0 & .y < 0)))

# duplicates
V10_sign |>
  dplyr::add_count(uniprot) |>
  dplyr::filter(n>1) |>
  dplyr::distinct()

# ## 101 - 6 = 95

# STEP 3: Top upregulated and downregulated genes
V10_sign %>%
  dplyr::filter(estimate_v1...v0_Durban > 0.95 & estimate_v1...v0_IAVI > 0.95)  # Top upregulated
V10_sign %>%
  dplyr::filter(estimate_v1...v0_Durban < -0.95 & estimate_v1...v0_IAVI < -0.95)  # Top downregulated


# STEP 4: Assign DEP_status based on log2FoldChange
V10_sign <- V10_sign %>%
  dplyr::mutate(DEP_status = ifelse(estimate_v1...v0_IAVI > 0, "UP", "DOWN"))

# STEP 5: Convert UniProt IDs to gene symbols using biomaRt
conversion_results_v10 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = V10_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# STEP 6: Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
V10_sign <- merge(V10_sign, conversion_results_v10, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

## top upregulated v10
V10_sign |> 
  dplyr::filter(estimate_v1...v0_Durban > 0.9  & estimate_v1...v0_IAVI > 0.9)

## top downregulated v10
V10_sign |> 
  dplyr::filter(estimate_v1...v0_Durban < -0.9  & estimate_v1...v0_IAVI < -0.9)

## number of proteins
V10_sign |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

## NO of down-regulated 
V10_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

## NO of up-regulated proteins
V10_sign |> 
  dplyr::filter(DEP_status == "UP") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

# STEP 7: Prepare data for visualization
longer_dfV10 <- V10_sign %>% 
  tidyr::pivot_longer(
    cols = estimate_v1...v0_IAVI:df_v1...v0_Durban,
    names_to = c(".value", "set", "cohort"),
    names_sep = "_"
  ) %>% 
  dplyr::mutate(
    ymax = estimate + qt(0.975, df = df) * SE,
    ymin = estimate - qt(0.975, df = df) * SE,
    ypdt = ymax * ymin,  # Identify intervals that do not include zero
    gene_dups = duplicated(Gene.Names)
  )

# STEP 8: Filter out genes with confidence intervals that include zero and remove non-replicated genes
longer_dfV10 <- longer_dfV10 %>%
  dplyr::filter(ypdt > 0) %>% # ypdt below zero gets all the confidence intervals that include a zero
  dplyr::filter(uniprot %in% unique(uniprot[duplicated(uniprot)])) # Keep only proteins appearing in both cohorts

# STEP 9: Remove Hemoglobin and immunoglobulin gegnes
longer_df_dupV10 <- longer_dfV10 %>%
  dplyr::filter(!grepl("^HB", Gene.Names) & !grepl("^IG", Gene.Names))  # Remove genes starting with 'HB' or 'IG'

# STEP 10: Retain only genes with strong effects in at least one cohort (|estimate| > 0.94) estimates > 0.94 or < -0.94
longer_df_dupV10 <- longer_df_dupV10 %>% 
  dplyr::group_by(Gene.Names) %>%
  dplyr::filter(any(estimate > 0.94 | estimate < -0.94)) %>%  # Keep only those where estimate > 1 in at least one cohort
  dplyr::ungroup()

# Step 11: Keep only genes appearing in both cohorts
longer_df_dup1V10 <- longer_df_dupV10 %>%
  dplyr::group_by(Gene.Names) %>%
  dplyr::filter(dplyr::n_distinct(cohort) > 1) %>%
  dplyr::ungroup()

# Step 12: Plotting the results
pdf("V10DEPS.pdf", width = 6, height = 3)
longer_df_dup1V10 |>  
  dplyr::mutate(name = forcats::fct_reorder(Gene.Names, estimate)) |>
  ggplot(aes(x=name, y = estimate, color = cohort , shape = exp))+
  geom_pointrange(aes(ymin = ymin, ymax = ymax),position=position_dodge(width=0.4), size=0.4) +
  geom_hline(yintercept = 0, col="gray",size =0.5, linetype = "dashed") + 
  scale_y_continuous(breaks=seq(-7,7,1), limits = c(-7,7))+
  scale_color_brewer(palette="Set1")+
  coord_flip()+
  theme_bw()+
  ggtitle("DEPs V1-V0") +theme(axis.line = element_line(color='black'),
                               plot.background = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.grid.major = element_blank()) 
dev.off() 


################################################################################
#### DIFFERENTIAL ANALYSIS
################################################################################
#B. Visit 2 - 0 Analysis

# STEP 1: Select columns of interest for V1-V0
V20 <- linear_df %>%
  dplyr::select(uniprot, PG.Genes, q.global, p.global, exp,
                estimate_v2...v0_IAVI, estimate_v2...v0_Durban,
                SE_v2...v0_IAVI, SE_v2...v0_Durban,
                df_v2...v0_IAVI, df_v2...v0_Durban)


# STEP 2: Filter for significant p-values and estimates in the same direction for both cohorts
V20_sign <- V20 %>%
  dplyr::filter(p.global < 0.005 & q.global < 0.005) %>%
  dplyr::filter(purrr::map2_lgl(estimate_v2...v0_Durban, estimate_v2...v0_IAVI, ~ (.x > 0 & .y > 0) | (.x < 0 & .y < 0)))

# duplicates
V20_sign |>
  dplyr::add_count(uniprot) |>
  dplyr::filter(n>1) |>
  dplyr::distinct()

# 182 - 14 = 168

# STEP 3: Top upregulated and downregulated genes
V20_sign %>%
  dplyr::filter(estimate_v2...v0_Durban > 0.95 & estimate_v2...v0_IAVI > 0.95)  # Top upregulated
V20_sign %>%
  dplyr::filter(estimate_v2...v0_Durban < -0.95 & estimate_v2...v0_IAVI < -0.95)  # Top downregulated


# STEP 4: Assign DEP_status based on log2FoldChange
V20_sign <- V20_sign %>%
  dplyr::mutate(DEP_status = ifelse(estimate_v2...v0_IAVI > 0, "UP", "DOWN"))

# STEP 5: Convert UniProt IDs to gene symbols using biomaRt

## Set up the biomaRt connection to Ensembl (human dataset) and check atrributes to use
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
biomaRt::listAttributes(mart = ensembl)

conversion_results_v20 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = V20_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# STEP 6: Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
V20_sign <- merge(V20_sign, conversion_results_v20, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

## top upregulated v20
V20_sign |> 
  dplyr::filter(estimate_v2...v0_Durban > 0.9  & estimate_v2...v0_IAVI > 0.9)

## top downregulated v20
V20_sign |> 
  dplyr::filter(estimate_v2...v0_Durban < -0.9  & estimate_v2...v0_IAVI < -0.9)

## number of proteins
V20_sign |> 
     dplyr::distinct(uniprot, .keep_all = TRUE) |>
     nrow()

## NO of down-regulated 
V20_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

## NO of up-regulated proteins
V20_sign |> 
  dplyr::filter(DEP_status == "UP") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()
# STEP 7: Prepare data for visualization
longer_dfV20 <- V20_sign %>% 
  tidyr::pivot_longer(
    cols = estimate_v2...v0_IAVI:df_v2...v0_Durban,
    names_to = c(".value", "set", "cohort"),
    names_sep = "_"
    ) %>% 
  dplyr::mutate(
    ymax = estimate + qt(0.975, df = df) * SE,
    ymin = estimate - qt(0.975, df = df) * SE,
    ypdt = ymax * ymin,  # Identify intervals that do not include zero
    gene_dups = duplicated(Gene.Names)
  )

# STEP 8: Filter out genes with confidence intervals that include zero and remove non-replicated genes
longer_dfV20 <- longer_dfV20 %>%
  dplyr::filter(ypdt > 0) %>% # ypdt below zero gets all the confidence intervals that include a zero
  dplyr::filter(uniprot %in% unique(uniprot[duplicated(uniprot)])) # Keep only proteins appearing in both cohorts

# STEP 9: Remove Hemoglobin and immunoglobulin gegnes
longer_df_dupV20 <- longer_dfV20 %>%
  dplyr::filter(!grepl("^HB", Gene.Names) & !grepl("^IG", Gene.Names))  # Remove genes starting with 'HB' or 'IG'

# STEP 10: Retain only genes with strong effects in at least one cohort (|estimate| > 0.94) estimates > 0.94 or < -0.94
longer_df_dupV20 <- longer_df_dupV20 %>% 
  dplyr::group_by(Gene.Names) %>%
  dplyr::filter(any(estimate > 0.94 | estimate < -0.94)) %>%  # Keep only those where estimate > 1 in at least one cohort
  dplyr::ungroup()

# Step 11: Keep only genes appearing in both cohorts
longer_df_dup1V20 <- longer_df_dupV20 %>%
  dplyr::group_by(Gene.Names) %>%
  dplyr::filter(dplyr::n_distinct(cohort) > 1) %>%
  dplyr::ungroup()

# Step 12: Plotting the results
pdf("V20DEPS.pdf", width = 6, height = 3)
longer_df_dup1V20 |>  
  dplyr::mutate(name = forcats::fct_reorder(Gene.Names, estimate)) |>
  ggplot(aes(x=name, y = estimate, color = cohort , shape = exp))+
  geom_pointrange(aes(ymin = ymin, ymax = ymax),position=position_dodge(width=0.4), size=0.4) +
  geom_hline(yintercept = 0, col="gray",size =0.5, linetype = "dashed") + 
  scale_y_continuous(breaks=seq(-7,7,1), limits = c(-7,7))+
  scale_color_brewer(palette="Set1")+
  coord_flip()+
  theme_bw()+
  ggtitle("DEPs V2-V0") +theme(axis.line = element_line(color='black'),
                                                         plot.background = element_blank(),
                                                         panel.grid.minor = element_blank(),
                                                         panel.grid.major = element_blank()) 
dev.off() 

################################################################################
#### DIFFERENTIAL ANALYSIS
################################################################################
#C. Visit 2 - 1 Analysis

# STEP 1: Select columns of interest for V2-V1
V21 <- linear_df %>%
  dplyr::select(uniprot, PG.Genes, q.global, p.global, exp,
                estimate_v2...v1_IAVI, estimate_v2...v1_Durban,
                SE_v2...v1_IAVI, SE_v2...v1_Durban,
                df_v2...v1_IAVI, df_v2...v1_Durban)


# STEP 2: Filter for significant p-values and estimates in the same direction for both cohorts
V21_sign <- V21 %>%
  dplyr::filter(p.global < 0.005 & q.global < 0.005) %>%
  dplyr::filter(purrr::map2_lgl(estimate_v2...v1_Durban, estimate_v2...v1_IAVI, ~ (.x > 0 & .y > 0) | (.x < 0 & .y < 0)))

# duplicates
V21_sign |>
  dplyr::add_count(uniprot) |>
  dplyr::filter(n>1) |>
  dplyr::distinct()

# number of unique DEPS at V21
V21_sign |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

# 162 - 13 = 149

# STEP 3: Top upregulated and downregulated genes
V21_sign %>%
  dplyr::filter(estimate_v2...v1_Durban > 0.95 & estimate_v2...v1_IAVI > 0.95)  # Top upregulated
V21_sign %>%
  dplyr::filter(estimate_v2...v1_Durban < -0.95 & estimate_v2...v1_IAVI < -0.95)  # Top downregulated


# STEP 4: Assign DEP_status based on log2FoldChange
V21_sign <- V21_sign %>%
  dplyr::mutate(DEP_status = ifelse(estimate_v2...v1_IAVI > 0, "UP", "DOWN"))

# STEP 5: Convert UniProt IDs to gene symbols using biomaRt

conversion_results_v21 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = V21_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# STEP 6: Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
V21_sign <- merge(V21_sign, conversion_results_v21, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

## top upregulated v21
V21_sign |> 
  dplyr::filter(estimate_v2...v1_Durban > 0.9  & estimate_v2...v1_IAVI > 0.9)

## top downregulated v20
V21_sign |> 
  dplyr::filter(estimate_v2...v1_Durban < -0.9  & estimate_v2...v1_IAVI < -0.9)

## NO of down-regulated 
V21_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

## NO of up-regulated proteins
V21_sign |> 
  dplyr::filter(DEP_status == "UP") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

# STEP 7: Prepare data for visualization
longer_dfV21 <- V21_sign %>% 
  tidyr::pivot_longer(
    cols = estimate_v2...v1_IAVI:df_v2...v1_Durban,
    names_to = c(".value", "set", "cohort"),
    names_sep = "_"
  ) %>% 
  dplyr::mutate(
    ymax = estimate + qt(0.975, df = df) * SE,
    ymin = estimate - qt(0.975, df = df) * SE,
    ypdt = ymax * ymin,  # Identify intervals that do not include zero
    gene_dups = duplicated(Gene.Names)
  )

# STEP 8: Filter out genes with confidence intervals that include zero and remove non-replicated genes
longer_dfV21 <- longer_dfV21 %>%
  dplyr::filter(ypdt > 0) %>% # ypdt below zero gets all the confidence intervals that include a zero
  dplyr::filter(uniprot %in% unique(uniprot[duplicated(uniprot)])) # Keep only proteins appearing in both cohorts

# STEP 9: Remove Hemoglobin and immunoglobulin gegnes
longer_df_dupV21 <- longer_dfV21 %>%
  dplyr::filter(!grepl("^HB", Gene.Names) & !grepl("^IG", Gene.Names))  # Remove genes starting with 'HB' or 'IG'

# STEP 10: Retain only genes with strong effects in at least one cohort (|estimate| > 0.94) estimates > 0.94 or < -0.94
longer_df_dupV21 <- longer_df_dupV21 %>% 
  dplyr::group_by(Gene.Names) %>%
  dplyr::filter(any(estimate > 0.94 | estimate < -0.94)) %>%  # Keep only those where estimate > 1 in at least one cohort
  dplyr::ungroup()

# Step 11: Keep only genes appearing in both cohorts
longer_df_dup1V21 <- longer_df_dupV21 %>%
  dplyr::group_by(Gene.Names) %>%
  dplyr::filter(dplyr::n_distinct(cohort) > 1) %>%
  dplyr::ungroup()

# Step 12: Plotting the results
pdf("V21DEPS.pdf", width = 6, height = 3)
longer_df_dup1V21 |>  
  dplyr::mutate(name = forcats::fct_reorder(Gene.Names, estimate)) |>
  ggplot(aes(x=name, y = estimate, color = cohort , shape = exp))+
  geom_pointrange(aes(ymin = ymin, ymax = ymax),position=position_dodge(width=0.4), size=0.4) +
  geom_hline(yintercept = 0, col="gray",size =0.5, linetype = "dashed") + 
  scale_y_continuous(breaks=seq(-7,7,1), limits = c(-7,7))+
  scale_color_brewer(palette="Set1")+
  coord_flip()+
  theme_bw()+
  ggtitle("DEPs V2-V1") +theme(axis.line = element_line(color='black'),
                               plot.background = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.grid.major = element_blank()) 
dev.off() 


################################################################################
#### PROTEIN INFORMATION
################################################################################

# A. V10
# Get protein function information for V10
V10_dep_function <- UniprotR::GetProteinFunction(V10_sign$uniprot)
V10_dep_function$ID <- rownames(V10_dep_function)

# Join V10_sign_0 with protein function data
V10_dep_function <- V10_sign %>%
  dplyr::inner_join(V10_dep_function, by = c("uniprot" = "ID"))

# Get subcellular location information for V10
V10_dep_loc <- UniprotR::GetSubcellular_location(V10_sign$uniprot)
V10_dep_loc$ID <- rownames(V10_dep_loc)

# Join protein function data with subcellular location data
V10_all_info <- V10_dep_function %>%
  dplyr::inner_join(V10_dep_loc, by = c("uniprot" = "ID"))

# Write the combined data to a file
data.table::fwrite(V10_all_info, "hiv_V10_deps_all_info.txt", sep = "\t")

# Retrieve known HIV-1 interaction proteins
hiv_interaction <- readr::read_csv("HIV-1_Interactions.csv")
hiv_interaction_V10_deps <- V10_sign  %>% 
  dplyr::inner_join(hiv_interaction, by = c("Gene.Names" = "Human_GeneSymbol"))

data.table::fwrite(hiv_interaction_V10_deps, "hiv_interaction_V10_deps.txt", sep = "\t")

#######-------------------------------------------------------------------------

# B. V20
# Get protein function information for V20
V20_dep_function <- UniprotR::GetProteinFunction(V20_sign$uniprot)
V20_dep_function$ID <- rownames(V20_dep_function)

# Join V20_sign_0 with protein function data
V20_dep_function <- V20_sign %>%
  dplyr::inner_join(V20_dep_function, by = c("uniprot" = "ID"))

# Get subcellular location information for V20
V20_dep_loc <- UniprotR::GetSubcellular_location(V20_sign$uniprot)
V20_dep_loc$ID <- rownames(V20_dep_loc)

# Join protein function data with subcellular location data
V20_all_info <- V20_dep_function %>%
  dplyr::inner_join(V20_dep_loc, by = c("uniprot" = "ID"))

# Write the combined data to a file
data.table::fwrite(V20_all_info, "hiv_V20_deps_all_info.txt", sep = "\t")

# Retrieve known HIV-1 interaction proteins
hiv_interaction <- readr::read_csv("HIV-1_Interactions.csv")
hiv_interaction_V20_deps <- V20_sign  %>% 
  dplyr::inner_join(hiv_interaction, by = c("Gene.Names" = "Human_GeneSymbol"))

data.table::fwrite(hiv_interaction_V20_deps, "hiv_interaction_V20_deps.txt", sep = "\t")


#######-------------------------------------------------------------------------

# C. V21
# Get protein function information for V21
V21_dep_function <- UniprotR::GetProteinFunction(V21_sign$uniprot)
V21_dep_function$ID <- rownames(V21_dep_function)

# Join V21_sign with protein function data
V21_dep_function <- V21_sign %>%
  dplyr::inner_join(V21_dep_function, by = c("uniprot" = "ID"))

# Get subcellular location information for V21
V21_dep_loc <- UniprotR::GetSubcellular_location(V21_sign$uniprot)
V21_dep_loc$ID <- rownames(V21_dep_loc)

# Join protein function data with subcellular location data
V21_all_info <- V21_dep_function %>%
  dplyr::inner_join(V21_dep_loc, by = c("uniprot" = "ID"))

# Write the combined data to a file
data.table::fwrite(V21_all_info, "hiv_V21_deps_all_info.txt", sep = "\t")

# Retrieve known HIV-1 interaction proteins
hiv_interaction <- readr::read_csv("HIV-1_Interactions.csv")
hiv_interaction_V21_deps <- V21_sign  %>% 
  dplyr::inner_join(hiv_interaction, by = c("Gene.Names" = "Human_GeneSymbol"))

data.table::fwrite(hiv_interaction_V21_deps, "hiv_interaction_V21_deps.txt", sep = "\t")



####################################
#OVER REPRESENTATION ANALYSIS
###################################

# Background data
backgrd_new <- readr::read_tsv("background_20221107.txt")

## NO of up-regulated proteins
V10_up <- V10_sign |> 
  dplyr::filter(DEP_status == "UP") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) 

# Step 1: Perform GO Enrichment Analysis
ego <- clusterProfiler::enrichGO(
  gene = V10_up$uniprot,
  universe = backgrd_new$Protein,
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Display enriched GO terms
head(ego)

# Extract significant results
cluster_profiler_v10up <- data.frame(ego@result)
padj0.2_v10up <- cluster_profiler_v10up %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)

# Simplify GO terms to remove redundancy
remove_v10up_redu <- clusterProfiler::simplify(ego)
cluster_p_v10_up <- as.data.frame(remove_v10up_redu@result)

####################################
# Preparations for Circos Plot
###################################

# Extract first gene name from Gene.Names column
V10_up <- V10_up %>%
  dplyr::mutate(Gene.Names = stringr::word(Gene.Names, 1, sep = "\\ ")) %>%
  dplyr::distinct(Gene.Names)

# Create a simplified data frame with gene descriptions
simplified_v10up <- cluster_p_v10_up %>%
  dplyr::select(geneID, Description)

# Initialize word count in gene names data frame
gene_names_v10up <- data.frame(keywords = V10_up$Gene.Names)
gene_names_v10up$word_count <- 0

# Populate word count for each process using a loop
for (i in seq_len(nrow(simplified_v10up))) {
  gene_names_v10up$word_count <- purrr::map_int(
    paste0("\\b", gene_names_v10up$keywords, "\\b"),
    ~ stringr::str_count(simplified_v10up$geneID[i], .x) %>% sum
  )
  colnames(gene_names_v10up)[i + 1] <- simplified_v10up$Description[i]
}

# Preview the gene names data frame
head(gene_names_v10up)

# Sum rows and columns for quantification
v10up_quant <- gene_names_v10up %>%
  tibble::column_to_rownames(var = "keywords") 

rsum10up <- rowSums(v10up_quant)
colSums(v10up_quant)
v10up_quant$rowsum <- rowSums(v10up_quant)
v10up_quant1<- gene_names_v10up %>%
  janitor::adorn_totals("row")

# Summarize by rowsum
v10up_quant %>%
  dplyr::group_by(rowsum) %>%
  dplyr::summarise(n = dplyr::n())


# Select relevant GO terms for further analysis
v10up_longmatrix2 <- data.frame(t(gene_names_v10up))
#v10up_longmatrix2[is.na(gene_names_v10up)] <- "igg"

colnames(v10up_longmatrix2) <- v10up_longmatrix2[1,] # convert the first row to column names
v10up_longmatrix2 <- v10up_longmatrix2[-1, ]   # remove the converted row names
#v21d_longmatrix2 <- cbind(rownames(v21d_longmatrix2), data.frame(v21d_longmatrix2, row.names=NULL)) # convert the row names to column

v10up_longmatrix2 <- readr::type_convert(v10up_longmatrix2)

# convert the df as a matrix
df.matrix_v10up <- as.matrix(v10up_longmatrix2, rownames.force = NA)


###############
## Hierarchical Clustering and Dendrogram to find/select the GO-TERMS

#optimal number of clusters
hier_v10up <- factoextra::fviz_nbclust(v10up_longmatrix2, kmeans, method = "s", k.max = 8)+
  ggplot2::ggtitle("optimal number of clusters for hierarchical clustering")
hier_v10up

##Basic dendogram

distMatrix_v10up <- dist(df.matrix_v10up, method = "euclidean")
groups_v10up <- hclust(distMatrix_v10up, method = "complete")
#groups_1$labels <- general_eugl$id


plot(groups_v10up, cex=0.8, hang=-1)
rect.hclust(groups_v10up, k=2)

##lad_id<- labels(groups_1)
dend_v10up <- as.dendrogram(groups_v10up)


pdf("v10up_dendogram.pdf", width = 12, height = 10)
par(mar=c(12, 4, 1, 4))
dend_v10up  %>% 
  dendextend::set("labels_cex", 0.7) %>% # Change size
  dendextend::set("branches_k_color", k = 6) %>% 
  plot(main = "BP enriched terms for v10_up") # plot
dev.off()

### select terms based on dendogram above
gene_names_v10up2 <- gene_names_v10up %>%
  dplyr::select(keywords, `response to stress`, `innate immune response`,
         `regulation of vesicle-mediated transport`, `adaptive immune response`,
         `transport`,`protein maturation`)

####################################
# Convert Matrix to Long Format
###################################

# Transpose the matrix and clean up column names
v10up_longmatrix2 <- as.data.frame(t(gene_names_v10up2))
colnames(v10up_longmatrix2) <- v10up_longmatrix2[1, ]
v10up_longmatrix2 <- v10up_longmatrix2[-1, ]

# Convert the data frame to a matrix
df.matrix <- as.matrix(readr::type_convert(v10up_longmatrix2))


####################################
# Circular Plots with Protein Atlas Data
###################################

# Read the protein atlas data
atlas <- readr::read_tsv("proteinatlas.tsv")
df.genes <- V10_up %>% 
  dplyr::distinct(Gene.Names)

# Join with atlas data
genes_atlas <- df.genes %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

# Prepare for circos plot
atlas_only_10up_enriched <- gene_names_v10up %>%
  dplyr::inner_join(atlas, by = c("keywords" = "Gene"))

# Handle missing data
#missing_v10_atlas <- gene_names_v10up[!gene_names_v10up$keywords %in% atlas_only_10up_enriched$keywords, ]

# Process and classify secretome locations
df9 <- atlas_only_10up_enriched %>%  dplyr::select("keywords",
                                                   "Secretome location",  "Subcellular location") 

df9$"Secretome location"[df9$"Secretome location" == "Secreted to blood"] <- "Secreted"
df9$"Secretome location"[is.na(df9$"Secretome location")] <- "Leakage"
df9$"Secretome location"[df9$"Secretome location" != "Secreted"] <- "Leakage"

df9 <- df9 %>%
  dplyr::select(keywords, `Secretome location`)

# Prepare final data for circular plot
#df12 <- rbind(df9, data.frame(keywords = stringr::word(missing_v10_atlas$keywords, 1), `Secretome location` = c("Secreted", "Secreted", "Leakage", "Leakage", "Leakage", "Secreted", "Secreted")))
df12 <- df9[order(df9$`Secretome location`), ]

df12 %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

df121 <- df12 %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

union(rownames(df.matrix), colnames(df.matrix))


####################################
# Circular Plot for V10 Upregulation
###################################

# Set up circos parameters
circlize::circos.clear()

##set the gaps between sectors — 
#circlize::circos.par(gap.degree=c(rep(2,nrow(df.matrix)-1),20,rep(2,ncol(df.matrix)-1),20))

# Define the order of sectors and color them
order1= union(rownames(df.matrix), colnames(df.matrix))

grid.col <- c("aquamarine4", "cadetblue4", "darkolivegreen4", "deepskyblue4","green", "navy", rep("orange", 23), rep("dimgrey", 28))
names(grid.col) <- union(rownames(df.matrix), df12$keywords)

# Draw the chord diagram
circlize::chordDiagram(df.matrix)

pdf("circos_plotV10up2.pdf", width = 8, height = 10)
circlize::chordDiagram(df.matrix, grid.col=grid.col,order = order1, annotationTrack = "grid", preAllocateTracks = 1)

# Add labels and axis
circlize::circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim <- circlize::get.cell.meta.data("xlim")
  ylim <- circlize::get.cell.meta.data("ylim")
  sector.name <- circlize::get.cell.meta.data("sector.index")
  
  circlize::circos.text(mean(xlim), ylim[1] + 1.5, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
}, bg.border = NA)

dev.off()

circlize::circos.clear()
dev.off()


## -----------------------------------------------------------------------------

## NO of down-regulated proteins
V10_down <- V10_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(Gene.Names, .keep_all = TRUE) 

V10_sign |> 
    dplyr::filter(DEP_status == "UP") |> 
    dplyr::distinct(Gene.Names, .keep_all = TRUE) |>
    nrow()

# Step 1: Perform GO Enrichment Analysis
ego_v10d <- clusterProfiler::enrichGO(
  gene = V10_down$uniprot,
  universe = backgrd_new$Protein,
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Display enriched GO terms
head(ego_v10d)

# Extract significant results
cluster_profiler_v10d <- data.frame(ego_v10d@result)
padj0.2_v10d <- cluster_profiler_v10d %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)

# Simplify GO terms to remove redundancy
remove_v10d_redu <- clusterProfiler::simplify(ego_v10d)
cluster_p_v10_d <- as.data.frame(remove_v10d_redu@result)

####################################
# Preparations for Circos Plot
###################################

# Extract first gene name from Gene.Names column
V10_d <- V10_down %>%
  dplyr::mutate(Gene.Names = stringr::word(Gene.Names, 1, sep = "\\ ")) %>%
  dplyr::distinct(Gene.Names)

# Create a simplified data frame with gene descriptions
simplified_v10d <- padj0.2_v10d %>%
  dplyr::select(geneID, Description)

# Initialize word count in gene names data frame
gene_names_v10d <- data.frame(keywords = V10_d$Gene.Names)
gene_names_v10d$word_count <- 0

# Populate word count for each process using a loop
for (i in seq_len(nrow(simplified_v10d))) {
  gene_names_v10d$word_count <- purrr::map_int(
    paste0("\\b", gene_names_v10d$keywords, "\\b"),
    ~ stringr::str_count(simplified_v10d$geneID[i], .x) %>% sum
  )
  colnames(gene_names_v10d)[i + 1] <- simplified_v10d$Description[i]
}

# Preview the gene names data frame
head(gene_names_v10d)

gene_names_v10d <- gene_names_v10d %>%
  dplyr::mutate(keywords = tidyr::replace_na(keywords, "IGKV2-29"))
#gene_names_v10d$keywords %>% tidyr::replace_na("IGKV2-29")

# Sum rows and columns for quantification
v10d_quant <- gene_names_v10d %>%
  tibble::column_to_rownames(var = "keywords") 

rsum10d <- rowSums(v10d_quant)
colSums(v10d_quant)
v10d_quant$rowsum <- rowSums(v10d_quant)
v10d_quant1<- gene_names_v10d %>%
  janitor::adorn_totals("row")

# Summarize by rowsum
v10d_quant %>%
  dplyr::group_by(rowsum) %>%
  dplyr::summarise(n = dplyr::n())


# Select relevant GO terms for further analysis
v10d_longmatrix2 <- data.frame(t(gene_names_v10d))
#v10up_longmatrix2[is.na(gene_names_v10up)] <- "igg"

colnames(v10d_longmatrix2) <- v10d_longmatrix2[1,] # convert the first row to column names
v10d_longmatrix2 <- v10d_longmatrix2[-1, ]   # remove the converted row names
#v21d_longmatrix2 <- cbind(rownames(v21d_longmatrix2), data.frame(v21d_longmatrix2, row.names=NULL)) # convert the row names to column


v10d_longmatrix2 <- readr::type_convert(v10d_longmatrix2)
colSums(v10d_longmatrix2)

# convert the df as a matrix
df.matrix_v10d <- as.matrix(v10d_longmatrix2, rownames.force = NA)


###############
## Hierarchical Clustering and Dendrogram to find/select the GO-TERMS

#optimal number of clusters
hier_v10d <- factoextra::fviz_nbclust(v10d_longmatrix2, kmeans, method = "s", k.max = 4)+
  ggplot2::ggtitle("optimal number of clusters for hierarchical clustering")
hier_v10d

##Basic dendogram

distMatrix_v10d <- dist(df.matrix_v10d, method = "euclidean")
groups_v10d <- hclust(distMatrix_v10d, method = "complete")
#groups_1$labels <- general_eugl$id


plot(groups_v10d, cex=0.8, hang=-1)
rect.hclust(groups_v10d, k=2)

##lad_id<- labels(groups_1)
dend_v10d <- as.dendrogram(groups_v10d)


pdf("v10down_dendogram.pdf", width = 12, height = 10)
par(mar=c(12, 4, 1, 4))
dend_v10d  %>% 
  dendextend::set("labels_cex", 0.7) %>% # Change size
  dendextend::set("branches_k_color", k = 2) %>% 
  plot(main = "BP enriched terms for v10_down") # plot
dev.off()

### select terms based on dendogram above
gene_names_v10d2 <- gene_names_v10d %>%
  dplyr::select(keywords, `establishment of endothelial barrier`, `regulation of peptidase activity`,
                `negative regulation of hydrolase activity`)

####################################
# Convert Matrix to Long Format
###################################

# Transpose the matrix and clean up column names
v10d_longmatrix2 <- as.data.frame(t(gene_names_v10d2))
colnames(v10d_longmatrix2) <- v10d_longmatrix2[1, ]
v10d_longmatrix2 <- v10d_longmatrix2[-1, ]

# Convert the data frame to a matrix
df.matrix_d <- as.matrix(readr::type_convert(v10d_longmatrix2))


####################################
# Circular Plots with Protein Atlas Data
###################################

# Read the protein atlas data
atlas <- readr::read_tsv("proteinatlas.tsv")
df.genes_d <- V10_d %>% 
  dplyr::distinct(Gene.Names)

df.genes_d <- df.genes_d %>%
  dplyr::mutate(keywords = tidyr::replace_na(Gene.Names, "IGKV2-29"))

# Join with atlas data
genes_atlas_d <- df.genes_d %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

# Prepare for circos plot
atlas_only_10d_enriched <- gene_names_v10d %>%
  dplyr::inner_join(atlas, by = c("keywords" = "Gene"))

# Handle missing data
#missing_v10_atlas <- gene_names_v10up[!gene_names_v10up$keywords %in% atlas_only_10up_enriched$keywords, ]

# Process and classify secretome locations
df9d <- atlas_only_10d_enriched %>%  dplyr::select("keywords",
                                                   "Secretome location",  "Subcellular location") 

df9d$"Secretome location"[df9d$"Secretome location" == "Secreted to blood"] <- "Secreted"
df9d$"Secretome location"[is.na(df9d$"Secretome location")] <- "Leakage"
df9d$"Secretome location"[df9d$"Secretome location" != "Secreted"] <- "Leakage"

df9d <- df9d %>%
  dplyr::select(keywords, `Secretome location`)

# Prepare final data for circular plot
#df12 <- rbind(df9, data.frame(keywords = stringr::word(missing_v10_atlas$keywords, 1), `Secretome location` = c("Secreted", "Secreted", "Leakage", "Leakage", "Leakage", "Secreted", "Secreted")))
df12d <- df9d[order(df9d$`Secretome location`), ]

df12d %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

df121d <- df12d %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

union(rownames(df.matrix_d), colnames(df.matrix_d))


####################################
# Circular Plot for V10 Upregulation
###################################

# Set up circos parameters
circlize::circos.clear()

##set the gaps between sectors — 
#circlize::circos.par(gap.degree=c(rep(2,nrow(df.matrix)-1),20,rep(2,ncol(df.matrix)-1),20))

# Define the order of sectors and color them
order1d= union(rownames(df.matrix_d), colnames(df.matrix_d))

grid.cold <- c("aquamarine4", "cadetblue4", "darkolivegreen4", rep("orange", 28), rep("dimgrey", 18))
names(grid.cold) <- union(rownames(df.matrix_d), colnames(df.matrix_d))

# Draw the chord diagram
circlize::chordDiagram(df.matrix_d)

pdf("circos_plotV10d2.pdf", width = 8, height = 10)
circlize::chordDiagram(df.matrix_d, grid.col=grid.cold,order = order1d, annotationTrack = "grid", preAllocateTracks = 1)

# Add labels and axis
circlize::circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim <- circlize::get.cell.meta.data("xlim")
  ylim <- circlize::get.cell.meta.data("ylim")
  sector.name <- circlize::get.cell.meta.data("sector.index")
  
  circlize::circos.text(mean(xlim), ylim[1] + 1.5, sector.name,
                        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)
}, bg.border = NA)

dev.off()

circlize::circos.clear()
dev.off()

################################################
## TISSUE DAMAGE
################################################


## GSEA with the tissue damage library
Tissue_damage_library <- readxl::read_excel("Tissue_damage_library.xlsx")

## Convert UniProt IDs to gene symbols using biomaRt
conversion_results_v10damage <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = V10$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

## Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
V10_damage <- merge(V10, conversion_results_v10damage, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)


## feature 1: numeric vector
geneList_v10 = V10_damage[,6]

## feature 2: named vector
names(geneList_v10) = as.character(V10_damage[,12])

## feature 3: decreasing order
geneList_v10 = sort(geneList_v10, decreasing = TRUE)

# Tissue damage enrichment

Tissue_gsea_v10 <- clusterProfiler::GSEA(geneList=geneList_v10,
                        TERM2GENE = Tissue_damage_library,
                        exponent = 1,
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps = 1e-10,
                        pvalueCutoff = 0.20,
                        pAdjustMethod = "none",
                        verbose = TRUE,
                        seed = TRUE)


cluster_profiler_GSEAv10t <- data.frame (Tissue_gsea_v10@result)
data.table::fwrite(cluster_profiler_GSEAv10t, "data/cluster_profiler_GSEAv10t.txt", sep = "\t")

pdf("tissue_damage_v10.pdf", width = 12, height = 10)
clusterProfiler::dotplot(Tissue_gsea_v10, showCategory = 10, title = "Visit 1 - Visit 0 associated Tissue damage signatures", split=".sign") + facet_grid(.~.sign)
dev.off()



##### secreted vs leakage 
# Read the protein atlas data
atlas <- readr::read_tsv("proteinatlas.tsv")
genes.all <- linear_df %>% 
  dplyr::distinct(uniprot)

# Join with atlas data
genes_atlas_all <- genes.all %>%
  dplyr::inner_join(atlas, by = c("uniprot" = "Uniprot"))

# Process and classify secretome locations
dfall <- genes_atlas_all %>%  dplyr::select("uniprot","Gene",
                                                   "Secretome location",  "Subcellular location") 

dfall$"Secretome location"[dfall$"Secretome location" == "Secreted to blood"] <- "Secreted"
dfall$"Secretome location"[is.na(dfall$"Secretome location")] <- "Leakage"
#df9d$"Secretome location"[df9d$"Secretome location" != "Secreted"] <- "Leakage"


dfall %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())
