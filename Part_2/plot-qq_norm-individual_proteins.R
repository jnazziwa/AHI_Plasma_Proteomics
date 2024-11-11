# -----------------------------------------------------------------------------#
# QQ plot to visually check the distribution of individual protein profiles
# -----------------------------------------------------------------------------#
# initiated on 2022-01-28
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())
# Please note some functions are called directly using '::' or ':::'
library(tidyverse)

#----- File names --------------------------------------------------------------

fn <- list(
  i = list(                               #  input
    nt01 = "../data/s1_neat_Spectronaut.v01.RData",
    dp01 = "../data/s1_depl_Spectronaut.v01.RData",
    interim = "../data/interim/Spectronaut_after_trimming.RData"
  ),
  o = list(                               #  output
    raw = "../results/figures/qq_norm-individual_proteins-raw.pdf"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))

#----- MAIN --------------------------------------------------------------------

# Load and combine proteomic data
e <- lapply(
  setNames(c("nt01", "dp01"), c("Neat", "Depleted")), 
  \(.x) {
    load(fn$i[[.x]], e <- new.env())
    e
  }
)
spct0 <- imap_dfr(e, ~mutate(.x$spct, prep_type = .y)) |> 
  mutate(protein_type_id = paste(protein_id, prep_type, sep = "-")) |> 
  separate(sampleid, c("subjid", "visit_nr"), sep = "-", remove = FALSE)
rm(e)

# Load trimmed data
load(fn$i$interim)


# IDs after trimming
after_trim <- distinct(spct, sampleid, protein_type_id, kept = "yes")

# kept = if each value is included after trimming of not
spct1 <- left_join(spct0, after_trim, by = c("sampleid", "protein_type_id")) |> 
  mutate(kept = replace_na(kept, "no")) |> 
  group_by(protein_type_id) |> 
  arrange(LogIntensities) |> 
  mutate(qn = qnorm(ppoints(n()))) |> 
  ungroup()


# Plotting ----------------------------------------------------------------

# too many plots so split
grouped_ids <- distinct(spct0, prep_type, protein_type_id) |> 
  arrange(prep_type) |> 
  pull(protein_type_id) |> 
  {\(x) split(x, ceiling(seq_along(x)/16))}()

rg <- range(spct1$LogIntensities)

pdf(fn$o$raw, height = 10, paper = "a4", pointsize = 10)

for(ea in grouped_ids) {
  p <- spct1 |> 
    filter(protein_type_id %in% ea) |> 
    ggplot() +
    geom_point(aes(qn, LogIntensities, color = kept), na.rm = T, size = 1) +
    scale_color_manual(guide = "none", 
                       values = c("no" = "red3", "yes" = "black")) +
    scale_y_continuous(limits = rg) +
    labs(x = "Theoretical", y = "Observed") +
    facet_wrap("protein_type_id", ncol = 4)
  print(p)
}
dev.off()
