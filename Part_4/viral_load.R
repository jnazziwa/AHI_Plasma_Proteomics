########
#VIRAL LOAD
###########


####################################
# VIRAL LOAD
####################################

# Read VL data for V20 and filter for significant proteins (p < 0.005)

VL_new <- xlsx::read.xlsx("lmer_VLvisit.xlsx", 1)
VL_newv10 <- VL_new |> dplyr::filter(visit_revpairwise == "v1 - v0")
VL_newv10_sign <- VL_newv10 |> dplyr::filter(p.global < 0.005  & p.value < 0.05)

VL_newv20 <- VL_new |> dplyr::filter(visit_revpairwise == "v2 - v0")
VL_newv20_sign <- VL_newv20 |> dplyr::filter(p.global < 0.005  & p.value < 0.05)

VL_newv21 <- VL_new |> dplyr::filter(visit_revpairwise == "v2 - v1")
VL_newv21_sign <- VL_newv21 |> dplyr::filter(p.global < 0.005  & p.value < 0.05)

VL_new_sign <- VL_new |> dplyr::filter(p.global < 0.005 & p.value < 0.05)
VL_new_sign$uniprot <- gsub("\\..*","",VL_new_sign$Protein)


####################################
# PIRATE Plots for VL V10
####################################

# Remove double accessions
VL_newv10_sign$uniprot <- gsub("\\..*","",VL_newv10_sign$Protein)

#  Convert UniProt IDs to gene symbols using biomaRt
conversion_results_VLv10 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = VL_newv10_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerVL_v10 <- merge(VL_newv10_sign, conversion_results_VLv10, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

# Create forest plot for VL V10
pdf("VL_V10_forestplot.pdf", width = 6, height = 2)
lmerVL_v10 %>%
  mutate(name = fct_reorder(Gene.Names, estimate)) %>%
  ggplot(aes(x = name, y = estimate, color = exp)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.4), size = 0.4) +
  geom_hline(yintercept = 0, col = "gray", size = 0.5, linetype = "dashed") +
  scale_y_continuous(breaks = seq(-7, 7, 1), limits = c(-7, 7)) +
  coord_flip() +
  labs(
    title = "Plot of Proteins Associated with VL at V1-V0 (p < 0.005)",
    x = "Protein", y = "log2FC"
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
dev.off()


####################################
# pIRATE Plots for VL V20
####################################

# Remove double accessions
VL_newv20_sign$uniprot <- gsub("\\..*","",VL_newv20_sign$Protein)

#  Convert UniProt IDs to gene symbols using biomaRt
conversion_results_VLv20 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = VL_newv20_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerVL_v20 <- merge(VL_newv20_sign, conversion_results_VLv20, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

# Create forest plot for VL V10
pdf("VL_V20_forestplot.pdf", width = 6, height = 2)
lmerVL_v20 %>%
  mutate(name = fct_reorder(Gene.Names, estimate)) %>%
  ggplot(aes(x = name, y = estimate, color = exp)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.4), size = 0.4) +
  geom_hline(yintercept = 0, col = "gray", size = 0.5, linetype = "dashed") +
  scale_y_continuous(breaks = seq(-7, 7, 1), limits = c(-7, 7)) +
  coord_flip() +
  labs(
    title = "Plot of Proteins Associated with VL at V2-V0 (p < 0.005)",
    x = "Protein", y = "log2FC"
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
dev.off()


####################################
# pIRATE Plots for VL V21
####################################

# Remove double accessions
VL_newv21_sign$uniprot <- gsub("\\..*","",VL_newv21_sign$Protein)

#  Convert UniProt IDs to gene symbols using biomaRt
conversion_results_VLv21 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = VL_newv21_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerVL_v21 <- merge(VL_newv21_sign, conversion_results_VLv21, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

# Create forest plot for VL V10
pdf("VL_V21_forestplot.pdf", width = 6, height = 2)
lmerVL_v21 %>%
  mutate(name = fct_reorder(Gene.Names, estimate)) %>%
  ggplot(aes(x = name, y = estimate, color = exp)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.4), size = 0.4) +
  geom_hline(yintercept = 0, col = "gray", size = 0.5, linetype = "dashed") +
  scale_y_continuous(breaks = seq(-7, 7, 1), limits = c(-7, 7)) +
  coord_flip() +
  labs(
    title = "Plot of Proteins Associated with VL at V2-V1 (p < 0.005)",
    x = "Protein", y = "log2FC"
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
dev.off()


####################################
# OVERREPRESENTATATION for VL 
####################################

#  Convert UniProt IDs to gene symbols using biomaRt
conversion_results_VL <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = VL_new_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerVL_sign <- merge(VL_new_sign, conversion_results_VL, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

lmerVL_sign$Gene.Names[lmerVL_sign$uniprot == 'P00736'] <- "C1R"
lmerVL_sign$Gene.Names[lmerVL_sign$uniprot == 'P02746'] <- "C1QB"

# Step 1: Perform GO Enrichment Analysis
ego_VL <- clusterProfiler::enrichGO(
  gene = lmerVL_sign$uniprot,
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
head(ego_VL)

# Extract significant results
cluster_profiler_VL <- data.frame(ego_VL@result)
padj0.2_VL <- cluster_profiler_VL %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)




################################################################################
#### PROTEIN INFORMATION
################################################################################

# A. VL
# Get protein function information for VL
lmerVL_sign <- lmerVL_sign %>% dplyr::distinct(uniprot, .keep_all = TRUE)
VL_function <- UniprotR::GetProteinFunction(lmerVL_sign$uniprot)
VL_function$ID <- rownames(VL_function)


# Join ESTIMATE with protein function data
VL_function <- lmerVL_sign %>%
  dplyr::inner_join(VL_function, by = c("uniprot" = "ID"))

# Get subcellular location information for V20
VL_loc <- UniprotR::GetSubcellular_location(lmerVL_sign$uniprot)
VL_loc$ID <- rownames(VL_loc)

# Join protein function data with subcellular location data
VL_all_info <- VL_function %>%
  dplyr::inner_join(VL_loc, by = c("uniprot" = "ID"))

# Write the combined data to a file
data.table::fwrite(VL_all_info, "VL_all_info.txt", sep = "\t")

# Retrieve known HIV-1 interaction proteins
hiv_interaction <- readr::read_csv("HIV-1_Interactions.csv")
hiv_interaction_VL <- lmerVL_sign  %>% 
  dplyr::inner_join(hiv_interaction, by = c("Gene.Names" = "Human_GeneSymbol"))

data.table::fwrite(hiv_interaction_VL, "hiv_interaction_VL.txt", sep = "\t")

#######-------------------------------------------------------------------------

# Join with atlas data
genes_atlas_vl <- lmerVL_sign %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

prt_info_VL <- genes_atlas_vl  %>%  dplyr::select("Gene.Names",
                                                    "Secretome location",  "Subcellular location", "Uniprot", "Protein class", "Biological process", "RNA tissue specific nTPM", "Gene description")
data.table::fwrite(prt_info_vl, "prt_info_vl.txt", sep = "\t")

