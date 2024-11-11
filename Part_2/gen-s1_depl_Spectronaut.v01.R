# -----------------------------------------------------------------------------#
# 1. Read and parse the output of Spectronaut
# 2. Make a list of UniProt ID for conversion table
# -----------------------------------------------------------------------------#
# initiated on 2021-06-29
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#

rm(list = ls())

# Please note some functions are called directly using '::' or ':::'
library(dplyr)
library(tidyr)
library(readr)    # read_csv, cols, col_*
stopifnot(getRversion() >= "4.1.0")   # R base pipe

#----- File names --------------------------------------------------------------

fn <- list(
  i = list(                               #  input
    spct = "../data/raw_internal/2022-02-01/depleted_20220201.txt.gz"
  ),
  o = list(                               #  output
    dp01 = "../data/s1_depl_Spectronaut.v01.RData"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))
# Warn any difference in input file(s) using MD5 hash
if(!all(
  tools::md5sum(fn$i$spct) == "b722092421308a4ea789ada468f42fb8"
)) warning(
  "One input file or more is not exactly identical to the one used during ",
  "development of this script."
)

#----- MAIN --------------------------------------------------------------------

# Read the Spectronaut file
spct0 <- read_csv(
  fn$i$spct,
  col_types = cols(
    Protein = col_character(),
    PG.Genes = col_character(),
    PG.ProteinNames = col_character(),
    intensity = col_double(),
    sampleID = col_character()
  )
) |> 
  drop_na(intensity)

# likely log transformed already
print(summary(spct0$intensity))

# As the format of the input file changed, the table was standardized like
# previous versions in order to fit other tables

spct <- spct0 |> 
  separate(sampleID, c("subjid", "visit_nr"), sep = "_") |> 
  mutate(subjid = sprintf("%03d", as.numeric(subjid))) |> 
  unite("sampleid", c(subjid, visit_nr), sep = "-") |> 
  mutate(
    # syntactically valid names
    protein_id = make.names(Protein),
    # likely log transformed already
    LogIntensities = intensity
  )

# Proteins
proteins <- spct |> 
  distinct(id = protein_id, Protein, PG.Genes, PG.ProteinNames)

##  Save -------------------------------------------------------------------

# original longitudinal format with more information
spct <- spct |> 
  select(sampleid, protein_id, LogIntensities)

save(proteins, spct, file= fn$o$dp01)
