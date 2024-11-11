
################################################################################
#### DIFFERENTIAL ANALYSIS
################################################################################

#A. Visit 2 - 0 Analysis

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

# ## 182 - 14 = 168

# STEP 3: Top upregulated and downregulated genes
V20_sign %>%
  dplyr::filter(estimate_v2...v0_Durban > 0.95 & estimate_v2...v0_IAVI > 0.95)  # Top upregulated
V20_sign %>%
  dplyr::filter(estimate_v2...v0_Durban < -0.95 & estimate_v2...v0_IAVI < -0.95)  # Top downregulated


# STEP 4: Assign DEP_status based on log2FoldChange
V20_sign <- V20_sign %>%
  dplyr::mutate(DEP_status = ifelse(estimate_v2...v0_IAVI > 0, "UP", "DOWN"))

# STEP 5: Convert UniProt IDs to gene symbols using biomaRt
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
#### PROTEIN INFORMATION
################################################################################

# A. V20
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



####################################
#OVER REPRESENTATION ANALYSIS
###################################

# Background data
backgrd_new <- readr::read_tsv("background_20221107.txt")

## NO of up-regulated proteins
V20_up <- V20_sign |> 
  dplyr::filter(DEP_status == "UP") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) 

V20_up <- V20_up %>%
  dplyr::mutate(Gene.Names = tidyr::replace_na(Gene.Names, "AMY1"))

# Step 1: Perform GO Enrichment Analysis
ego_v20up <- clusterProfiler::enrichGO(
  gene = V20_up$uniprot,
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
head(ego_v20up)

# Extract significant results
cluster_profiler_v20up <- data.frame(ego_v20up@result)
padj0.2_v20up <- cluster_profiler_v20up %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)

# Simplify GO terms to remove redundancy
remove_v20up_redu <- clusterProfiler::simplify(ego_v20up)
cluster_p_v20_up <- as.data.frame(remove_v20up_redu@result)

####################################
# Preparations for Circos Plot
###################################

# Extract first gene name from Gene.Names column
V20_up <- V20_up %>%
  dplyr::mutate(Gene.Names = stringr::word(Gene.Names, 1, sep = "\\ ")) %>%
  dplyr::distinct(Gene.Names)

# Create a simplified data frame with gene descriptions
simplified_v20up <- cluster_p_v20_up %>%
  dplyr::select(geneID, Description)

# Initialize word count in gene names data frame
gene_names_v20up <- data.frame(keywords = V20_up$Gene.Names)
gene_names_v20up$word_count <- 0

# Populate word count for each process using a loop
for (i in seq_len(nrow(simplified_v20up))) {
  gene_names_v20up$word_count <- purrr::map_int(
    paste0("\\b", gene_names_v20up$keywords, "\\b"),
    ~ stringr::str_count(simplified_v20up$geneID[i], .x) %>% sum
  )
  colnames(gene_names_v20up)[i + 1] <- simplified_v20up$Description[i]
}

# Preview the gene names data frame
head(gene_names_v20up)

# Sum rows and columns for quantification
v20up_quant <- gene_names_v20up %>%
  tibble::column_to_rownames(var = "keywords") 

rsum20up <- rowSums(v20up_quant)
colSums(v20up_quant)
v20up_quant$rowsum <- rowSums(v20up_quant)
v20up_quant1<- gene_names_v20up %>%
  janitor::adorn_totals("row")

# Summarize by rowsum
v20up_quant %>%
  dplyr::group_by(rowsum) %>%
  dplyr::summarise(n = dplyr::n())


# Select relevant GO terms for further analysis
v20up_longmatrix2 <- data.frame(t(gene_names_v20up))
#v20up_longmatrix2[is.na(gene_names_v20up)] <- "igg"

colnames(v20up_longmatrix2) <- v20up_longmatrix2[1,] # convert the first row to column names
v20up_longmatrix2 <- v20up_longmatrix2[-1, ]   # remove the converted row names
#v21d_longmatrix2 <- cbind(rownames(v21d_longmatrix2), data.frame(v21d_longmatrix2, row.names=NULL)) # convert the row names to column

v20up_longmatrix2 <- readr::type_convert(v20up_longmatrix2)
colSums(v20up_longmatrix2)

# convert the df as a matrix
df.matrix_v20up <- as.matrix(v20up_longmatrix2, rownames.force = NA)


###############
## Hierarchical Clustering and Dendrogram to find/select the GO-TERMS

#optimal number of clusters
hier_v20up <- factoextra::fviz_nbclust(v20up_longmatrix2, kmeans, method = "s", k.max = 4)+
  ggplot2::ggtitle("optimal number of clusters for hierarchical clustering")
hier_v20up

##Basic dendogram

distMatrix_v20up <- dist(df.matrix_v20up, method = "euclidean")
groups_v20up <- hclust(distMatrix_v20up, method = "complete")
#groups_1$labels <- general_eugl$id


plot(groups_v20up, cex=0.8, hang=-1)
rect.hclust(groups_v20up, k=2)

##lad_id<- labels(groups_1)
dend_v20up <- as.dendrogram(groups_v20up)


pdf("v20up_dendogram.pdf", width = 12, height = 10)
par(mar=c(12, 4, 1, 4))
dend_v20up  %>% 
  dendextend::set("labels_cex", 0.7) %>% # Change size
  dendextend::set("branches_k_color", k = 6) %>% 
  plot(main = "BP enriched terms for v20_up") # plot
dev.off()

### select terms based on dendogram above
# gene_names_v20up2 <- gene_names_v20up %>%
#   dplyr::select(keywords, `response to stress`, `innate immune response`,
#                 `regulation of vesicle-mediated transport`, `adaptive immune response`,
#                 `transport`,`protein maturation`)

####################################
# Convert Matrix to Long Format
###################################

# Transpose the matrix and clean up column names
v20up_longmatrix2 <- as.data.frame(t(gene_names_v20up))
colnames(v20up_longmatrix2) <- v20up_longmatrix2[1, ]
v20up_longmatrix2 <- v20up_longmatrix2[-1, ]

# Convert the data frame to a matrix
df.matrix20 <- as.matrix(readr::type_convert(v20up_longmatrix2))


####################################
# Circular Plots with Protein Atlas Data
###################################

# Read the protein atlas data
atlas <- readr::read_tsv("proteinatlas.tsv")
df.genes20 <- V20_up %>% 
  dplyr::distinct(Gene.Names)

# Join with atlas data
genes_atlas20 <- df.genes20 %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

# Prepare for circos plot
atlas_only_20up_enriched <- gene_names_v20up %>%
  dplyr::inner_join(atlas, by = c("keywords" = "Gene"))

prt_info_V20UP <- atlas_only_20up_enriched  %>%  dplyr::select("keywords",
                                                         "Secretome location",  "Subcellular location", "Uniprot", "Protein class", "Biological process", "RNA tissue specific nTPM", "Gene description")
data.table::fwrite(prt_info_V20UP, "prt_info_V20UP.txt", sep = "\t")



# Handle missing data
#missing_v20_atlas <- gene_names_v20up[!gene_names_v20up$keywords %in% atlas_only_20up_enriched$keywords, ]

# Process and classify secretome locations
df920 <- atlas_only_20up_enriched %>%  dplyr::select("keywords",
                                                   "Secretome location",  "Subcellular location") 

df920$"Secretome location"[df920$"Secretome location" == "Secreted to blood"] <- "Secreted"
df920$"Secretome location"[is.na(df920$"Secretome location")] <- "Leakage"
df920$"Secretome location"[df920$"Secretome location" != "Secreted"] <- "Leakage"

df920 <- df920 %>%
  dplyr::select(keywords, `Secretome location`)

# Prepare final data for circular plot
#df12 <- rbind(df9, data.frame(keywords = stringr::word(missing_v10_atlas$keywords, 1), `Secretome location` = c("Secreted", "Secreted", "Leakage", "Leakage", "Leakage", "Secreted", "Secreted")))
df1220 <- df920[order(df920$`Secretome location`), ]

df1220 %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

df12120 <- df1220 %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

union(rownames(df.matrix20), colnames(df.matrix20))


####################################
# Circular Plot for V20 Upregulation
###################################

# Set up circos parameters
circlize::circos.clear()

##set the gaps between sectors — 
#circlize::circos.par(gap.degree=c(rep(2,nrow(df.matrix)-1),20,rep(2,ncol(df.matrix)-1),20))

# Define the order of sectors and color them
order120= union(rownames(df.matrix20), colnames(df.matrix20))

grid.col20 <- c("aquamarine4", "cadetblue4", "darkolivegreen4", "deepskyblue4","green", "navy", rep("orange", 58), rep("dimgrey", 18))
names(grid.col20) <-  union(rownames(df.matrix20), df1220$keywords)

# Draw the chord diagram
circlize::chordDiagram(df.matrix20)

pdf("circos_plotV20up2.pdf", width = 8, height = 10)
circlize::chordDiagram(df.matrix20, grid.col=grid.col20 ,order = order120, annotationTrack = "grid", preAllocateTracks = 1)

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
V20_down <- V20_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(Gene.Names, .keep_all = TRUE) 

V20_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(Gene.Names, .keep_all = TRUE) |>
  nrow()

# Step 1: Perform GO Enrichment Analysis
ego_v20d <- clusterProfiler::enrichGO(
  gene = V20_down$uniprot,
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
head(ego_v20d)

# Extract significant results
cluster_profiler_v20d <- data.frame(ego_v20d@result)
padj0.2_v20d <- cluster_profiler_v20d %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)

# Simplify GO terms to remove redundancy
remove_v20d_redu <- clusterProfiler::simplify(ego_v20d)
cluster_p_v20_d <- as.data.frame(remove_v20d_redu@result)

####################################
# Preparations for Circos Plot
###################################

# Extract first gene name from Gene.Names column
V20_d <- V20_down %>%
  dplyr::mutate(Gene.Names = stringr::word(Gene.Names, 1, sep = "\\ ")) %>%
  dplyr::distinct(Gene.Names)

# Create a simplified data frame with gene descriptions
simplified_v20d <- padj0.2_v20d %>%
  dplyr::select(geneID, Description)

# Initialize word count in gene names data frame
gene_names_v20d <- data.frame(keywords = V20_d$Gene.Names)
gene_names_v20d$word_count <- 0

# Populate word count for each process using a loop
for (i in seq_len(nrow(simplified_v20d))) {
  gene_names_v20d$word_count <- purrr::map_int(
    paste0("\\b", gene_names_v20d$keywords, "\\b"),
    ~ stringr::str_count(simplified_v20d$geneID[i], .x) %>% sum
  )
  colnames(gene_names_v20d)[i + 1] <- simplified_v20d$Description[i]
}

# Preview the gene names data frame
head(gene_names_v20d)

gene_names_v20d <- gene_names_v20d %>%
  dplyr::mutate(keywords = tidyr::replace_na(keywords, "IGKV2-29"))
#gene_names_v20d$keywords %>% tidyr::replace_na("IGKV2-29")

# Sum rows and columns for quantification
v20d_quant <- gene_names_v20d %>%
  tibble::column_to_rownames(var = "keywords") 

rsum20d <- rowSums(v20d_quant)
colSums(v20d_quant)
v20d_quant$rowsum <- rowSums(v20d_quant)
v20d_quant1<- gene_names_v20d %>%
  janitor::adorn_totals("row")

# Summarize by rowsum
v20d_quant %>%
  dplyr::group_by(rowsum) %>%
  dplyr::summarise(n = dplyr::n())


# Select relevant GO terms for further analysis
v20d_longmatrix2 <- data.frame(t(gene_names_v20d))
#v20up_longmatrix2[is.na(gene_names_v20up)] <- "igg"

colnames(v20d_longmatrix2) <- v20d_longmatrix2[2,] # convert the first row to column names
v20d_longmatrix2 <- v20d_longmatrix2[-1, ]   # remove the converted row names
#v21d_longmatrix2 <- cbind(rownames(v21d_longmatrix2), data.frame(v21d_longmatrix2, row.names=NULL)) # convert the row names to column


v20d_longmatrix2 <- readr::type_convert(v20d_longmatrix2)
colSums(v20d_longmatrix2)

# convert the df as a matrix
df.matrix_v20d <- as.matrix(v20d_longmatrix2, rownames.force = NA)


###############
## Hierarchical Clustering and Dendrogram to find/select the GO-TERMS

#optimal number of clusters
hier_v20d <- factoextra::fviz_nbclust(v20d_longmatrix2, kmeans, method = "s", k.max = 4)+
  ggplot2::ggtitle("optimal number of clusters for hierarchical clustering")
hier_v20d

##Basic dendogram

distMatrix_v20d <- dist(df.matrix_v20d, method = "euclidean")
groups_v20d <- hclust(distMatrix_v20d, method = "complete")
#groups_1$labels <- general_eugl$id


plot(groups_v20d, cex=0.8, hang=-1)
rect.hclust(groups_v20d, k=7)

##lad_id<- labels(groups_1)
dend_v20d <- as.dendrogram(groups_v20d)


pdf("v20down_dendogram.pdf", width = 12, height = 10)
par(mar=c(12, 4, 1, 4))
dend_v20d  %>% 
  dendextend::set("labels_cex", 0.7) %>% # Change size
  dendextend::set("branches_k_color", k = 8) %>% 
  plot(main = "BP enriched terms for v20_down") # plot
dev.off()

### select terms based on dendogram above
gene_names_v20d2 <- gene_names_v20d %>%
  dplyr::select(keywords,`regulation of peptidase activity`,`inflammatory response`,`response to stress`, 
                `positive regulation of wound healing`, `regulation of immune response`, `tumor necrosis factor superfamily cytokine production`,
                `neutrophil chemotaxis`)

####################################
# Convert Matrix to Long Format
###################################

# Transpose the matrix and clean up column names
v20d_longmatrix2 <- as.data.frame(t(gene_names_v20d2))
colnames(v20d_longmatrix2) <- v20d_longmatrix2[1, ]
v20d_longmatrix2 <- v20d_longmatrix2[-1, ]

# Convert the data frame to a matrix
df.matrix_20d <- as.matrix(readr::type_convert(v20d_longmatrix2))


####################################
# Circular Plots with Protein Atlas Data
###################################

# Read the protein atlas data
atlas <- readr::read_tsv("proteinatlas.tsv")
df.genes_20d <- V20_d %>% 
  dplyr::distinct(Gene.Names)

df.genes_20d <- df.genes_20d %>%
  dplyr::mutate(keywords = tidyr::replace_na(Gene.Names, "IGKV2-29"))

# Join with atlas data
genes_atlas_20d <- df.genes_20d %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

# Prepare for circos plot
atlas_only_20d_enriched <- gene_names_v20d %>%
  dplyr::inner_join(atlas, by = c("keywords" = "Gene"))

prt_info_V20d <- atlas_only_20d_enriched  %>%  dplyr::select("keywords",
                                                               "Secretome location",  "Subcellular location", "Uniprot", "Protein class", "Biological process", "RNA tissue specific nTPM", "Gene description")
data.table::fwrite(prt_info_V20d, "prt_info_V20d.txt", sep = "\t")


# Handle missing data
#missing_v20_atlas <- gene_names_v20up[!gene_names_v20up$keywords %in% atlas_only_20up_enriched$keywords, ]

# Process and classify secretome locations
df920d <- atlas_only_20d_enriched %>%  dplyr::select("keywords",
                                                   "Secretome location",  "Subcellular location") 

df920d$"Secretome location"[df920d$"Secretome location" == "Secreted to blood"] <- "Secreted"
df920d$"Secretome location"[is.na(df920d$"Secretome location")] <- "Leakage"
df920d$"Secretome location"[df920d$"Secretome location" != "Secreted"] <- "Leakage"

df920d <- df920d %>%
  dplyr::select(keywords, `Secretome location`)

# Prepare final data for circular plot
#df12 <- rbind(df9, data.frame(keywords = stringr::word(missing_v10_atlas$keywords, 1), `Secretome location` = c("Secreted", "Secreted", "Leakage", "Leakage", "Leakage", "Secreted", "Secreted")))
df1220d <- df920d[order(df920d$`Secretome location`), ]

df1220d %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

df12120d <- df1220d %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

union(rownames(df.matrix_20d), colnames(df.matrix_20d))


####################################
# Circular Plot for V20 Upregulation
###################################

# Set up circos parameters
circlize::circos.clear()

##set the gaps between sectors — 
#circlize::circos.par(gap.degree=c(rep(2,nrow(df.matrix)-1),20,rep(2,ncol(df.matrix)-1),20))

# Define the order of sectors and color them
order120d= union(rownames(df.matrix_20d), colnames(df.matrix_20d))

grid.col20d <- c("aquamarine4", "cadetblue4", "darkolivegreen4","deepskyblue4","green", "navy", "maroon", rep("orange", 48), rep("dimgrey", 44))
names(grid.col20d) <- union(rownames(df.matrix_20d), df1220d$keywords)

# Draw the chord diagram
circlize::chordDiagram(df.matrix_20d)

pdf("circos_plotV20d2.pdf", width = 8, height = 10)
circlize::chordDiagram(df.matrix_20d, grid.col=grid.col20d,order = order120d, annotationTrack = "grid", preAllocateTracks = 1)

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
conversion_results_v20damage <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = V20$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

## Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
V20_damage <- merge(V20, conversion_results_v20damage, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)


## feature 1: numeric vector
geneList_v20 = V20_damage[,6]

## feature 2: named vector
names(geneList_v20) = as.character(V20_damage[,12])

## feature 3: decreasing order
geneList_v20 = sort(geneList_v20, decreasing = TRUE)

# Tissue damage enrichment

Tissue_gsea_v20 <- clusterProfiler::GSEA(geneList=geneList_v20,
                                         TERM2GENE = Tissue_damage_library,
                                         exponent = 1,
                                         minGSSize = 10,
                                         maxGSSize = 500,
                                         eps = 1e-10,
                                         pvalueCutoff = 0.20,
                                         pAdjustMethod = "none",
                                         verbose = TRUE,
                                         seed = TRUE)


cluster_profiler_GSEAv20t <- data.frame (Tissue_gsea_v20@result)
data.table::fwrite(cluster_profiler_GSEAv20t, "data/cluster_profiler_GSEAv20t.txt", sep = "\t")

pdf("tissue_damage_v20.pdf", width = 12, height = 10)
clusterProfiler::dotplot(Tissue_gsea_v20, showCategory = 10, title = "Visit 2 - Visit 0 associated Tissue damage signatures", split=".sign") + facet_grid(.~.sign)
dev.off()