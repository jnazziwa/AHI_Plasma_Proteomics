# -----------------------------------------------------------------------------#
# Read the table and create RData file of clinical info
# -----------------------------------------------------------------------------#
# initiated on 2020-11-13
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#

rm(list = ls())

# Please note some functions are called directly using '::' or ':::'
library(dplyr)
#library(tidyr)
library(readr)    # read_csv, cols, col_*
stopifnot(getRversion() >= "4.1.0")   # R base pipe

#----- File names --------------------------------------------------------------

fn <- list(
  i = list(                               #  input
    pat = "../data/raw_internal/2020-10-23/Patient_Matched_sampleIDs_updated_AH_02042020 (003).xlsx",
    ars01 = "../data/raw_internal/2021-01-11/symptoms_to_Mun-Gwan_11012021.csv",  # ARS
    vl01 = "../data/s1-viral_load.v01.RData"      # viral load
  ),
  o = list(                               #  output
    c01 = "../data/s1-clinic.v01.RData"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(fn$lib, unlist(fn$i), dirname(unlist(fn$o)))))
# Warn any difference in input file(s) using MD5 hash
if(!all(
  tools::md5sum(fn$i$pat) == "ee15f59e4a93709d47d706ee7888144f",
  tools::md5sum(fn$i$ars01) == "97438139c144dbbc8f299080aed5613a",
  tools::md5sum(fn$i$vl01) == "95f47739ec9c07e17a2dd00dc6e705e7"
)) warning(
  "One input file or more is not exactly identical to the one used during ",
  "development of this script."
)

#----- MAIN --------------------------------------------------------------------


##  Read clinical info table -----------------------------------------------
#  Clinical data per patient
clinc0 <- readxl::read_xlsx(
  fn$i$pat, sheet = "demo_subtype_ars",
  col_types = c("text", "text", "text", "date", "date", rep("text", 13))
) 

#  read ARS data of IAVI samples
ars <- read_csv(
  fn$i$ars01,
  col_types = cols_only(
    pid = col_character(),
    ars = col_character()
  )
) |> 
  mutate(ars = factor(ars, levels = c("No", "Yes")))

stopifnot(
  anyDuplicated(ars$pid) == 0,    # unique pids
  all(ars$pid %in% clinc0$pid)    # no mismatch
)

clinc <- clinc0 |> 
  left_join(ars, by = "pid") |>    # add ARS

  ##  Make it tidy
  #  e.g `fatigue1` -> `Fatigue`
  rename_with(~tools::toTitleCase(sub("1$", "", .x))) |> 
  #  Yes/No variables
  mutate(across(
    Fever:Anorexia,     # symptoms
    ~ factor(.x, levels = c("0", "1"), labels = c("No", "Yes"))
  )) |> 
  rename(
    "subjid" = "Pid",  #  to make patient IDs identical to other tables.
    "ARS" = "Ars"
  ) |> 
  mutate(
    subjid = if_else(nchar(subjid) > 3,
                     subjid,
                     formatC(as.numeric(subjid), width = 3, flag = "0")),
    Site = tools::toTitleCase(Site),
    Sex = factor(Sex, levels= c("F", "M"), labels= c("Female", "Male")),
    Cohort = if_else(Site == "Durban", "Durban", "IAVI")
  )


##  Add rough estimate of VL set point -------------------------------------

#  Load viral load data
load(fn$i$vl01)

#  VL set point computation
vl_setpt <- left_join(vload, clinc, by= "subjid") |> 
  mutate(
    #  days from EDI
    day_fr_edi = ((visitdate - Edi) / 3600 / 24)      # sec -> day
  ) |> 
  group_by(subjid) |> 
  summarise(
    #  median of viral load, 3 weeks after EDI
    vl_setpt = median(vload[day_fr_edi > 21], na.rm= T),
    .groups = "drop"
  )

#  confirm no missing patient  
stopifnot(nrow(vl_setpt) == nrow(clinc))

#  update `clinc` by adding VL setpt
clinc <- clinc |> 
  right_join(vl_setpt, ., by= "subjid")


##  Save -------------------------------------------------------------------

save(clinc, file= fn$o$c01)
