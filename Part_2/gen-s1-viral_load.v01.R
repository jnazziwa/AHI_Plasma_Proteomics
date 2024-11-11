# -----------------------------------------------------------------------------#
# Read the table and create RData file of viral load
# -----------------------------------------------------------------------------#
# initiated on 2020-11-12
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#

rm(list = ls())

# Please note some functions are called directly using '::' or ':::'
library(dplyr)
stopifnot(getRversion() >= "4.1.0")   # R base pipe

#----- File names --------------------------------------------------------------

fn <- list(
  i = list(                               #  input
    pat = "../data/raw_internal/2020-10-23/Patient_Matched_sampleIDs_updated_AH_02042020 (003).xlsx"
  ),
  o = list(                               #  output
    vl01 = "../data/s1-viral_load.v01.RData"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))
# Warn any difference in input file(s) using MD5 hash
if(!all(
  tools::md5sum(fn$i$pat) == "ee15f59e4a93709d47d706ee7888144f"
)) warning(
  "One input file or more is not exactly identical to the one used during ",
  "development of this script."
)

#----- MAIN --------------------------------------------------------------------

#   viral load, CD4, and CD8 per sample (or visit)
vload0 <- readxl::read_xlsx(
  fn$i$pat,
  sheet = "vload_cd4_cd8",
  col_types = c("text", "numeric", "numeric", "numeric", "date", "text")
) 


##  Save -------------------------------------------------------------------

vload <- vload0 |> 
  mutate(
    Timepoint = factor(Timepoint, levels = c("0", "1", "2"))
  ) |> 

  #  to make patient IDs through all tables with patient info.
  rename("subjid" = "pid") |> 
  mutate(
    subjid = if_else(
      nchar(subjid) > 3,
      subjid,
      formatC(as.numeric(subjid), width = 3, flag = "0")   # e.g. "036" in Durban
    )
  )

save(vload, file= fn$o$vl)
