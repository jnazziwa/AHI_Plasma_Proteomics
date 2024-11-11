########
#ARS
###########


####################################
# ARS
####################################

# Read ARS data for V20 and filter for significant proteins (p < 0.005)

ARS_new <- xlsx::read.xlsx("lmer_ARSvisit.xlsx", 1)
ARS_newv10 <- ARS_new |> dplyr::filter(visit_revpairwise == "v1 - v0")
ARS_newv10_sign <- ARS_newv10 |> dplyr::filter(p.global < 0.005  & p.value < 0.05)

ARS_newv20 <- ARS_new |> dplyr::filter(visit_revpairwise == "v2 - v0")
ARS_newv20_sign <- ARS_newv20 |> dplyr::filter(p.global < 0.005  & p.value < 0.05)

ARS_new_sign <- ARS_new |> dplyr::filter(p.global < 0.005)
ARS_new_sign$uniprot <- gsub("\\..*","",ARS_new_sign$Protein)


####################################
# pIRATE Plots for ARS V10
####################################

# Remove double accessions
ARS_newv10_sign$uniprot <- gsub("\\..*","",ARS_newv10_sign$Protein)

#  Convert UniProt IDs to gene symbols using biomaRt
conversion_results_ARSv10 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = ARS_newv10_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerARS_v10 <- merge(ARS_newv10_sign, conversion_results_ARSv10, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

# Create forest plot for ARS V10
pdf("ARS_V10_forestplot.pdf", width = 6, height = 2)
lmerARS_v10 %>%
  mutate(name = fct_reorder(Gene.Names, estimate)) %>%
  ggplot(aes(x = name, y = estimate, color = exp)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.4), size = 0.4) +
  geom_hline(yintercept = 0, col = "gray", size = 0.5, linetype = "dashed") +
  scale_y_continuous(breaks = seq(-7, 7, 1), limits = c(-7, 7)) +
  coord_flip() +
  labs(
    title = "Plot of Proteins Associated with ARS at V1-V0 (p < 0.005)",
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
# pIRATE Plots for ARS V20
####################################

# Remove double accessions
ARS_newv20_sign$uniprot <- gsub("\\..*","",ARS_newv20_sign$Protein)

#  Convert UniProt IDs to gene symbols using biomaRt
conversion_results_ARSv20 <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = ARS_newv20_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerARS_v20 <- merge(ARS_newv20_sign, conversion_results_ARSv20, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

# Create forest plot for ARS V10
pdf("ARS_V20_forestplot.pdf", width = 6, height = 2)
lmerARS_v20 %>%
  mutate(name = fct_reorder(Gene.Names, estimate)) %>%
  ggplot(aes(x = name, y = estimate, color = exp)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.4), size = 0.4) +
  geom_hline(yintercept = 0, col = "gray", size = 0.5, linetype = "dashed") +
  scale_y_continuous(breaks = seq(-7, 7, 1), limits = c(-7, 7)) +
  coord_flip() +
  labs(
    title = "Plot of Proteins Associated with ARS at V2-V0 (p < 0.005)",
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
# OVERREPRESENTATATION for ARS 
####################################

#  Convert UniProt IDs to gene symbols using biomaRt
conversion_results_ARS <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = ARS_new_sign$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerARS_sign <- merge(ARS_new_sign, conversion_results_ARS, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)

lmerARS_sign$Gene.Names[lmerARS_sign$uniprot == 'P00736'] <- "C1R"
lmerARS_sign$Gene.Names[lmerARS_sign$uniprot == 'P02746'] <- "C1QB"

# Step 1: Perform GO Enrichment Analysis
ego_ARS <- clusterProfiler::enrichGO(
  gene = lmerARS_sign$uniprot,
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
head(ego_ARS)

# Extract significant results
cluster_profiler_ARS <- data.frame(ego_ARS@result)
padj0.2_ARS <- cluster_profiler_ARS %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)


####################################
# Preparations for Circos Plot
###################################


# Create a simplified data frame with gene descriptions
simplified_ARS <- padj0.2_ARS %>%
  dplyr::select(geneID, Description)

# Initialize word count in gene names data frame
gene_names_ARS <- data.frame(keywords = lmerARS_sign$Gene.Names)
gene_names_ARS$word_count <- 0

# Populate word count for each process using a loop
for (i in seq_len(nrow(simplified_ARS))) {
  gene_names_ARS$word_count <- purrr::map_int(
    paste0("\\b", gene_names_ARS$keywords, "\\b"),
    ~ stringr::str_count(simplified_ARS$geneID[i], .x) %>% sum
  )
  colnames(gene_names_ARS)[i + 1] <- simplified_ARS$Description[i]
}

# Preview the gene names data frame
head(gene_names_ARS)
gene_names_ARS <- gene_names_ARS %>% dplyr::distinct(keywords, .keep_all = TRUE)



gene_names_ARS <- gene_names_ARS %>%
  dplyr::mutate(keywords = tidyr::replace_na(keywords, "IGKV2-29"))
#gene_names_ARS$keywords %>% tidyr::replace_na("IGKV2-29")

# Sum rows and columns for quantification
ARS_quant <- gene_names_ARS %>%
  tibble::column_to_rownames(var = "keywords") 

rsumARS <- rowSums(ARS_quant)
colSums(ARS_quant)
ARS_quant$rowsum <- rowSums(ARS_quant)
ARS_quant1<- gene_names_ARS %>%
  janitor::adorn_totals("row")

# Summarize by rowsum
ARS_quant %>%
  dplyr::group_by(rowsum) %>%
  dplyr::summarise(n = dplyr::n())


# Select relevant GO terms for further analysis
ARS_longmatrix2 <- data.frame(t(gene_names_ARS))
#v21up_longmatrix2[is.na(gene_names_v21up)] <- "igg"

colnames(ARS_longmatrix2) <- ARS_longmatrix2[2,] # convert the first row to column names
ARS_longmatrix2 <- ARS_longmatrix2[-1, ]   # remove the converted row names
#ARS_longmatrix2 <- cbind(rownames(ARS_longmatrix2), data.frame(ARS_longmatrix2, row.names=NULL)) # convert the row names to column


ARS_longmatrix2 <- readr::type_convert(ARS_longmatrix2)
colSums(ARS_longmatrix2)

# convert the df as a matrix
df.matrix_ARS <- as.matrix(ARS_longmatrix2, rownames.force = NA)


###############
## Hierarchical Clustering and Dendrogram to find/select the GO-TERMS

#optimal number of clusters
hier_ARS <- factoextra::fviz_nbclust(ARS_longmatrix2, kmeans, method = "s", k.max = 4)+
  ggplot2::ggtitle("optimal number of clusters for hierarchical clustering")
hier_ARS

##Basic dendogram

distMatrix_ARS <- dist(df.matrix_ARS, method = "euclidean")
groups_ARS <- hclust(distMatrix_ARS, method = "complete")
#groups_1$labels <- general_eugl$id


plot(groups_ARS, cex=0.8, hang=-1)
rect.hclust(groups_ARS, k=4)

##lad_id<- labels(groups_1)
dend_ARS <- as.dendrogram(groups_ARS)


pdf("ARS_dendogram.pdf", width = 12, height = 10)
par(mar=c(12, 4, 1, 4))
dend_ARS  %>% 
  dendextend::set("labels_cex", 0.7) %>% # Change size
  dendextend::set("branches_k_color", k = 4) %>% 
  plot(main = "BP enriched terms for v21_down") # plot
dev.off()

################################################################################
#### PROTEIN INFORMATION
################################################################################

# A. ARS
# Get protein function information for ARS
lmerARS_sign <- lmerARS_sign %>% dplyr::distinct(uniprot, .keep_all = TRUE)
ARS_function <- UniprotR::GetProteinFunction(lmerARS_sign$uniprot)
ARS_function$ID <- rownames(ARS_function)


# Join V20_sign_0 with protein function data
ARS_function <- lmerARS_sign %>%
  dplyr::inner_join(ARS_function, by = c("uniprot" = "ID"))

# Get subcellular location information for V20
ARS_loc <- UniprotR::GetSubcellular_location(lmerARS_sign$uniprot)
ARS_loc$ID <- rownames(ARS_loc)

# Join protein function data with subcellular location data
ARS_all_info <- ARS_function %>%
  dplyr::inner_join(ARS_loc, by = c("uniprot" = "ID"))

# Write the combined data to a file
data.table::fwrite(ARS_all_info, "ARS_all_info.txt", sep = "\t")

# Retrieve known HIV-1 interaction proteins
hiv_interaction <- readr::read_csv("HIV-1_Interactions.csv")
hiv_interaction_ARS <- lmerARS_sign  %>% 
  dplyr::inner_join(hiv_interaction, by = c("Gene.Names" = "Human_GeneSymbol"))

data.table::fwrite(hiv_interaction_ARS, "hiv_interaction_ARS.txt", sep = "\t")

#######-------------------------------------------------------------------------



#### classifcation of patterns

meancluster8 <- readr::read_tsv("meancluster8.txt")
meancluster8$Protein <-  stringr::word(meancluster8$Protein,1,sep = "\\.")

# Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
lmerARS_gp <- merge(lmerARS_sign, meancluster8, by.x = "uniprot", by.y = "Protein", all.x = TRUE) 


# Join with atlas data
genes_atlas_ARS <- lmerARS_sign %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

prt_info_ars <- genes_atlas_ARS  %>%  dplyr::select("Gene.Names",
                                                             "Secretome location",  "Subcellular location", "Uniprot", "Protein class", "Biological process", "RNA tissue specific nTPM", "Gene description")
data.table::fwrite(prt_info_ars, "prt_info_ars.txt", sep = "\t")

