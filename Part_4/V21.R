
################################################################################
#### DIFFERENTIAL ANALYSIS
################################################################################

#A. Visit 2 - 1 Analysis

# STEP 1: Select columns of interest for V1-V2
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

# ## 162 - 13 = 149

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

## number of proteins
V21_sign |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) |>
  nrow()

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

# A. V21
# Get protein function information for V21
V21_dep_function <- UniprotR::GetProteinFunction(V21_sign$uniprot)
V21_dep_function$ID <- rownames(V21_dep_function)

# Join V21_sign_1 with protein function data
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

#######-------------------------------------------------------------------------



####################################
#OVER REPRESENTATION ANALYSIS
###################################

# Background data
backgrd_new <- readr::read_tsv("background_20221107.txt")

## NO of up-regulated proteins
V21_up <- V21_sign |> 
  dplyr::filter(DEP_status == "UP") |> 
  dplyr::distinct(uniprot, .keep_all = TRUE) 

V21_up <- V21_up %>%
  dplyr::mutate(Gene.Names = tidyr::replace_na(Gene.Names, "AMY1"))

# Step 1: Perform GO Enrichment Analysis
ego_v21up <- clusterProfiler::enrichGO(
  gene = V21_up$uniprot,
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
head(ego_v21up)

# Extract significant results
cluster_profiler_v21up <- data.frame(ego_v21up@result)
padj0.2_v21up <- cluster_profiler_v21up %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)

# Simplify GO terms to remove redundancy
remove_v21up_redu <- clusterProfiler::simplify(ego_v21up)
cluster_p_v21_up <- as.data.frame(remove_v21up_redu@result)

####################################
# Preparations for Circos Plot
###################################

# Extract first gene name from Gene.Names column
V21_up <- V21_up %>%
  dplyr::mutate(Gene.Names = stringr::word(Gene.Names, 1, sep = "\\ ")) %>%
  dplyr::distinct(Gene.Names)

# Create a simplified data frame with gene descriptions
simplified_v21up <- padj0.2_v21up %>%
  dplyr::select(geneID, Description)

# Initialize word count in gene names data frame
gene_names_v21up <- data.frame(keywords = V21_up$Gene.Names)
gene_names_v21up$word_count <- 0

# Populate word count for each process using a loop
for (i in seq_len(nrow(simplified_v21up))) {
  gene_names_v21up$word_count <- purrr::map_int(
    paste0("\\b", gene_names_v21up$keywords, "\\b"),
    ~ stringr::str_count(simplified_v21up$geneID[i], .x) %>% sum
  )
  colnames(gene_names_v21up)[i + 1] <- simplified_v21up$Description[i]
}

# Preview the gene names data frame
head(gene_names_v21up)

# Sum rows and columns for quantification
v21up_quant <- gene_names_v21up %>%
  tibble::column_to_rownames(var = "keywords") 

rsum21up <- rowSums(v21up_quant)
colSums(v21up_quant)
v21up_quant$rowsum <- rowSums(v21up_quant)
v21up_quant1<- gene_names_v21up %>%
  janitor::adorn_totals("row")

# Summarize by rowsum
v21up_quant %>%
  dplyr::group_by(rowsum) %>%
  dplyr::summarise(n = dplyr::n())


# Select relevant GO terms for further analysis
v21up_longmatrix2 <- data.frame(t(gene_names_v21up))
#v21up_longmatrix2[is.na(gene_names_v21up)] <- "igg"

colnames(v21up_longmatrix2) <- v21up_longmatrix2[1,] # convert the first row to column names
v21up_longmatrix2 <- v21up_longmatrix2[-1, ]   # remove the converted row names
#v21d_longmatrix2 <- cbind(rownames(v21d_longmatrix2), data.frame(v21d_longmatrix2, row.names=NULL)) # convert the row names to column

v21up_longmatrix2 <- readr::type_convert(v21up_longmatrix2)
colSums(v21up_longmatrix2)

# convert the df as a matrix
df.matrix_v21up <- as.matrix(v21up_longmatrix2, rownames.force = NA)


###############
## Hierarchical Clustering and Dendrogram to find/select the GO-TERMS

#optimal number of clusters
hier_v21up <- factoextra::fviz_nbclust(v21up_longmatrix2, kmeans, method = "s", k.max = 2)+
  ggplot2::ggtitle("optimal number of clusters for hierarchical clustering")
hier_v21up

##Basic dendogram

distMatrix_v21up <- dist(df.matrix_v21up, method = "euclidean")
groups_v21up <- hclust(distMatrix_v21up, method = "complete")
#groups_1$labels <- general_eugl$id


plot(groups_v21up, cex=0.8, hang=-1)
rect.hclust(groups_v21up, k=2)

##lad_id<- labels(groups_1)
dend_v21up <- as.dendrogram(groups_v21up)


pdf("v21up_dendogram.pdf", width = 12, height = 10)
par(mar=c(12, 4, 1, 4))
dend_v21up  %>% 
  dendextend::set("labels_cex", 0.7) %>% # Change size
  dendextend::set("branches_k_color", k = 3) %>% 
  plot(main = "BP enriched terms for v21_up") # plot
dev.off()

### select terms based on dendogram above
# gene_names_v21up2 <- gene_names_v21up %>%
#   dplyr::select(keywords, `response to stress`, `innate immune response`,
#                 `regulation of vesicle-mediated transport`, `adaptive immune response`,
#                 `transport`,`protein maturation`)

####################################
# Convert Matrix to Long Format
###################################

# Transpose the matrix and clean up column names
v21up_longmatrix2 <- as.data.frame(t(gene_names_v21up))
colnames(v21up_longmatrix2) <- v21up_longmatrix2[1, ]
v21up_longmatrix2 <- v21up_longmatrix2[-1, ]

# Convert the data frame to a matrix
df.matrix21 <- as.matrix(readr::type_convert(v21up_longmatrix2))


####################################
# Circular Plots with Protein Atlas Data
###################################

# Read the protein atlas data
atlas <- readr::read_tsv("proteinatlas.tsv")
df.genes21 <- V21_up %>% 
  dplyr::distinct(Gene.Names)

# Join with atlas data
genes_atlas21 <- df.genes21 %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

# Prepare for circos plot
atlas_only_21up_enriched <- gene_names_v21up %>%
  dplyr::inner_join(atlas, by = c("keywords" = "Gene"))

# Handle missing data
#missing_v21_atlas <- gene_names_v21up[!gene_names_v21up$keywords %in% atlas_only_21up_enriched$keywords, ]

# Process and classify secretome locations
df921 <- atlas_only_21up_enriched %>%  dplyr::select("keywords",
                                                     "Secretome location",  "Subcellular location") 

df921$"Secretome location"[df921$"Secretome location" == "Secreted to blood"] <- "Secreted"
df921$"Secretome location"[is.na(df921$"Secretome location")] <- "Leakage"
df921$"Secretome location"[df921$"Secretome location" != "Secreted"] <- "Leakage"

df921 <- df921 %>%
  dplyr::select(keywords, `Secretome location`)

# Prepare final data for circular plot
#df12 <- rbind(df9, data.frame(keywords = stringr::word(missing_v11_atlas$keywords, 1), `Secretome location` = c("Secreted", "Secreted", "Leakage", "Leakage", "Leakage", "Secreted", "Secreted")))
df1221 <- df921[order(df921$`Secretome location`), ]

df1221 %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

df12121 <- df1221 %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

union(rownames(df.matrix21), colnames(df.matrix21))


####################################
# Circular Plot for V21 Upregulation
###################################

# Set up circos parameters
circlize::circos.clear()

##set the gaps between sectors — 
#circlize::circos.par(gap.degree=c(rep(2,nrow(df.matrix)-1),20,rep(2,ncol(df.matrix)-1),20))

# Define the order of sectors and color them
order121= union(rownames(df.matrix21), colnames(df.matrix21))

grid.col21 <- c("aquamarine4", "cadetblue4", "darkolivegreen4", rep("orange", 49), rep("dimgrey", 13))
names(grid.col21) <-  union(rownames(df.matrix21), df1221$keywords)

# Draw the chord diagram
circlize::chordDiagram(df.matrix21)

pdf("circos_plotV21up2.pdf", width = 8, height = 10)
circlize::chordDiagram(df.matrix21, grid.col=grid.col21 ,order = order121, annotationTrack = "grid", preAllocateTracks = 1)

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
V21_down <- V21_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(Gene.Names, .keep_all = TRUE) 

V21_sign |> 
  dplyr::filter(DEP_status == "DOWN") |> 
  dplyr::distinct(Gene.Names, .keep_all = TRUE) |>
  nrow()

# Step 1: Perform GO Enrichment Analysis
ego_v21d <- clusterProfiler::enrichGO(
  gene = V21_down$uniprot,
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
head(ego_v21d)

# Extract significant results
cluster_profiler_v21d <- data.frame(ego_v21d@result)
padj0.2_v21d <- cluster_profiler_v21d %>% 
  dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)

# Simplify GO terms to remove redundancy
remove_v21d_redu <- clusterProfiler::simplify(ego_v21d)
cluster_p_v21_d <- as.data.frame(remove_v21d_redu@result)

####################################
# Preparations for Circos Plot
###################################

# Extract first gene name from Gene.Names column
V21_d <- V21_down %>%
  dplyr::mutate(Gene.Names = stringr::word(Gene.Names, 1, sep = "\\ ")) %>%
  dplyr::distinct(Gene.Names)

# Create a simplified data frame with gene descriptions
simplified_v21d <- padj0.2_v21d %>%
  dplyr::select(geneID, Description)

# Initialize word count in gene names data frame
gene_names_v21d <- data.frame(keywords = V21_d$Gene.Names)
gene_names_v21d$word_count <- 0

# Populate word count for each process using a loop
for (i in seq_len(nrow(simplified_v21d))) {
  gene_names_v21d$word_count <- purrr::map_int(
    paste0("\\b", gene_names_v21d$keywords, "\\b"),
    ~ stringr::str_count(simplified_v21d$geneID[i], .x) %>% sum
  )
  colnames(gene_names_v21d)[i + 1] <- simplified_v21d$Description[i]
}

# Preview the gene names data frame
head(gene_names_v21d)

gene_names_v21d <- gene_names_v21d %>%
  dplyr::mutate(keywords = tidyr::replace_na(keywords, "IGKV2-29"))
#gene_names_v21d$keywords %>% tidyr::replace_na("IGKV2-29")

# Sum rows and columns for quantification
v21d_quant <- gene_names_v21d %>%
  tibble::column_to_rownames(var = "keywords") 

rsum21d <- rowSums(v21d_quant)
colSums(v21d_quant)
v21d_quant$rowsum <- rowSums(v21d_quant)
v21d_quant1<- gene_names_v21d %>%
  janitor::adorn_totals("row")

# Summarize by rowsum
v21d_quant %>%
  dplyr::group_by(rowsum) %>%
  dplyr::summarise(n = dplyr::n())


# Select relevant GO terms for further analysis
v21d_longmatrix2 <- data.frame(t(gene_names_v21d))
#v21up_longmatrix2[is.na(gene_names_v21up)] <- "igg"

colnames(v21d_longmatrix2) <- v21d_longmatrix2[2,] # convert the first row to column names
v21d_longmatrix2 <- v21d_longmatrix2[-1, ]   # remove the converted row names
#v21d_longmatrix2 <- cbind(rownames(v21d_longmatrix2), data.frame(v21d_longmatrix2, row.names=NULL)) # convert the row names to column


v21d_longmatrix2 <- readr::type_convert(v21d_longmatrix2)
colSums(v21d_longmatrix2)

# convert the df as a matrix
df.matrix_v21d <- as.matrix(v21d_longmatrix2, rownames.force = NA)


###############
## Hierarchical Clustering and Dendrogram to find/select the GO-TERMS

#optimal number of clusters
hier_v21d <- factoextra::fviz_nbclust(v21d_longmatrix2, kmeans, method = "s", k.max = 4)+
  ggplot2::ggtitle("optimal number of clusters for hierarchical clustering")
hier_v21d

##Basic dendogram

distMatrix_v21d <- dist(df.matrix_v21d, method = "euclidean")
groups_v21d <- hclust(distMatrix_v21d, method = "complete")
#groups_1$labels <- general_eugl$id


plot(groups_v21d, cex=0.8, hang=-1)
rect.hclust(groups_v21d, k=7)

##lad_id<- labels(groups_1)
dend_v21d <- as.dendrogram(groups_v21d)


pdf("v21down_dendogram.pdf", width = 12, height = 10)
par(mar=c(12, 4, 1, 4))
dend_v21d  %>% 
  dendextend::set("labels_cex", 0.7) %>% # Change size
  dendextend::set("branches_k_color", k = 7) %>% 
  plot(main = "BP enriched terms for v21_down") # plot
dev.off()

### select terms based on dendogram above
gene_names_v21d2 <- gene_names_v21d %>%
  dplyr::select(keywords,`inflammatory response`,`response to stress`, `innate immune response`,
                `detoxification`, `response to biotic stimulus`, `complement activation`,`neutrophil chemotaxis`) 
                
                
####################################
# Convert Matrix to Long Format
###################################

# Transpose the matrix and clean up column names
v21d_longmatrix2 <- as.data.frame(t(gene_names_v21d2))
colnames(v21d_longmatrix2) <- v21d_longmatrix2[1, ]
v21d_longmatrix2 <- v21d_longmatrix2[-1, ]

# Convert the data frame to a matrix
df.matrix_21d <- as.matrix(readr::type_convert(v21d_longmatrix2))


####################################
# Circular Plots with Protein Atlas Data
###################################

# Read the protein atlas data
atlas <- readr::read_tsv("proteinatlas.tsv")
df.genes_21d <- V21_d %>% 
  dplyr::distinct(Gene.Names)

df.genes_21d <- df.genes_21d %>%
  dplyr::mutate(keywords = tidyr::replace_na(Gene.Names, "IGKV2-29"))

# Join with atlas data
genes_atlas_21d <- df.genes_21d %>%
  dplyr::inner_join(atlas, by = c("Gene.Names" = "Gene"))

# Prepare for circos plot
atlas_only_21d_enriched <- gene_names_v21d %>%
  dplyr::inner_join(atlas, by = c("keywords" = "Gene"))

# Handle missing data
#missing_v21_atlas <- gene_names_v21up[!gene_names_v21up$keywords %in% atlas_only_21up_enriched$keywords, ]

# Process and classify secretome locations
df921d <- atlas_only_21d_enriched %>%  dplyr::select("keywords",
                                                     "Secretome location",  "Subcellular location") 

df921d$"Secretome location"[df921d$"Secretome location" == "Secreted to blood"] <- "Secreted"
df921d$"Secretome location"[is.na(df921d$"Secretome location")] <- "Leakage"
df921d$"Secretome location"[df921d$"Secretome location" != "Secreted"] <- "Leakage"

df921d <- df921d %>%
  dplyr::select(keywords, `Secretome location`)

# Prepare final data for circular plot
#df12 <- rbind(df9, data.frame(keywords = stringr::word(missing_v11_atlas$keywords, 1), `Secretome location` = c("Secreted", "Secreted", "Leakage", "Leakage", "Leakage", "Secreted", "Secreted")))
df1221d <- df921d[order(df921d$`Secretome location`), ]

df1221d %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

df12121d <- df1221d %>% 
  dplyr::group_by(`Secretome location`) %>%
  dplyr::summarise(n=dplyr::n())

union(rownames(df.matrix_21d), colnames(df.matrix_21d))


####################################
# Circular Plot for V21 Upregulation
###################################

# Set up circos parameters
circlize::circos.clear()

##set the gaps between sectors — 
#circlize::circos.par(gap.degree=c(rep(2,nrow(df.matrix)-1),21,rep(2,ncol(df.matrix)-1),21))

# Define the order of sectors and color them
order121d= union(rownames(df.matrix_21d), colnames(df.matrix_21d))

grid.col21d <- c("aquamarine4", "cadetblue4", "darkolivegreen4","deepskyblue4", "navy", "maroon","darkred", rep("orange", 40), rep("dimgrey", 46))
names(grid.col21d) <- union(rownames(df.matrix_21d), df1221d$keywords)

# Draw the chord diagram
circlize::chordDiagram(df.matrix_21d)

pdf("circos_plotV21d2.pdf", width = 8, height = 10)
circlize::chordDiagram(df.matrix_21d, grid.col=grid.col21d,order = order121d, annotationTrack = "grid", preAllocateTracks = 1)

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
conversion_results_v21damage <- biomaRt::getBM(
  filters = "uniprotswissprot",
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  values = V21$uniprot,
  mart = ensembl
) %>%
  dplyr::distinct(uniprotswissprot, .keep_all = TRUE) #remove duplicates

## Merge gene symbols back into the data and rename HGNC symbol to Gene.Names
V21_damage <- merge(V21, conversion_results_v21damage, by.x = "uniprot", by.y = "uniprotswissprot", all.x = TRUE) %>% 
  dplyr::rename(Gene.Names = hgnc_symbol)


## feature 1: numeric vector
geneList_v21 = V21_damage[,6]

## feature 2: named vector
names(geneList_v21) = as.character(V21_damage[,12])

## feature 3: decreasing order
geneList_v21 = sort(geneList_v21, decreasing = TRUE)

# Tissue damage enrichment

Tissue_gsea_v21 <- clusterProfiler::GSEA(geneList=geneList_v21,
                                         TERM2GENE = Tissue_damage_library,
                                         exponent = 1,
                                         minGSSize = 10,
                                         maxGSSize = 500,
                                         eps = 1e-10,
                                         pvalueCutoff = 0.20,
                                         pAdjustMethod = "none",
                                         verbose = TRUE,
                                         seed = TRUE)


cluster_profiler_GSEAv21t <- data.frame (Tissue_gsea_v21@result)
data.table::fwrite(cluster_profiler_GSEAv21t, "data/cluster_profiler_GSEAv21t.txt", sep = "\t")

pdf("tissue_damage_v21.pdf", width = 12, height = 10)
clusterProfiler::dotplot(Tissue_gsea_v21, showCategory = 10, title = "Visit 2 - Visit 1 associated Tissue damage signatures", split=".sign") + facet_grid(.~.sign)
dev.off()
