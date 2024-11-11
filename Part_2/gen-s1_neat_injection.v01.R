# -----------------------------------------------------------------------------#
# 1. Read and parse the tjbles about injection
# 2. make `durban` clean according to info (meeting note : 2020-10-26)
# 3. make `injtn` tidy
#    a. split to patientID and visit number
#    b. attach `durban`
# 4. Make a list of UniProt ID for conversion table
# -----------------------------------------------------------------------------#
# initiated on 2022-01-31 split from `gen-s1_neat_Spectronaut.v01.R`
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

# Please note some functions are called directly using '::' or ':::'
library(dplyr)
stopifnot(getRversion() >= "4.1.0")   # R base pipe

#----- File names --------------------------------------------------------------

fn <- list(
  i = list(                               #  input
    injtn = "../data/raw_internal/2020-10-23/Patient_Matched_sampleIDs_updated_AH_02042020 (003).xlsx",
    durban = "../data/raw_internal/2020-10-23/05. Sample_infomation_Durban.xlsx"
  ),
  o = list(                               #  output
    injn = "../data/s1_neat_injection.v01.RData"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))
# Warn any difference in input file(s) using MD5 hash
if(!all(
  tools::md5sum(fn$i$injtn) == "ee15f59e4a93709d47d706ee7888144f",
  tools::md5sum(fn$i$durban) == "ef0dae5817624619efb058acf855d266"
)) warning(
  "One input file or more is not exactly identical to the one used during ",
  "development of this script."
)

#----- MAIN --------------------------------------------------------------------

#  Durban samples' info (e.g. VL-cps/ml, Days PI)
durban0 <- readxl::read_xlsx(fn$i$durban, "Samples from Durban")

#  Clean `durban`
durban <- durban0 |> 
  #  remove unnecessary columns
  select(-ID...4, -RandomisationID, -`Randomisation #`, -`Position In Box`,
         -`Added 14 days`, -`Time point`, -`Visit code`) |> 
  #  to date format assuming origin was set to '1900-01-01'
  mutate(
    `Date Collected` = as.POSIXct(`Date Collected` * 60 * 60 * 24, origin= "1900-01-01")
  ) |> 
  rename(
    ID = "ID...1", 
    run_id = "RunID",
    collect_date = "Date Collected"
  )

#  data per injection? including some others such as VL setpoint, Country
injtn0 <- readxl::read_xlsx(
  fn$i$injtn, "All_Injections_Info_First_Proje",
  col_types = "text"
) 

#  attach `durban`
injtn_d <- injtn0 |> 
  select(run_id = "File Name", PatientID) |> 
  left_join(durban, by= "run_id")
#  confirm the column `ID` was redundant
stopifnot(identical(injtn_d$PatientID[!is.na(injtn_d$ID)], 
                    injtn_d$ID[!is.na(injtn_d$ID)]))

##  make `injtn` tidy
injtn <- injtn_d |> 
  #  remove the redundant `ID` column
  select(-ID) |> 
  rename(sid_in_cohort = "PatientID") |> 
  #  extract 'subjid' = patient unique ID
  #          'visit_nr' = visit number
  tidyr::separate(sid_in_cohort, c("subjid", "visit_nr"), sep= "_") |> 
  mutate(
    Cohort = if_else(run_id %in% durban$run_id, "Durban", "IAVI"),
    subjid = stringr::str_extract(subjid, "[:digit:]+$"),
    # to integer
    visit_nr = case_when(
      visit_nr == "Screening" ~ "1",
      visit_nr == "B02" ~ "2",
      !is.na(VisitNumber) ~ format(VisitNumber),   # durban
      TRUE ~ visit_nr
    ) |> 
      as.integer()
  ) |> 
  mutate(sampleid = paste0(subjid, "-T", visit_nr), .before = run_id) |> 
  select(-VisitNumber)    # duplicated

#  confirm all `RunID`s in `injtn` are in `RunID` column of `durban0` tibble
durbanids_in_injtn_vl <- injtn$run_id[startsWith(injtn$run_id, "1910")]
stopifnot(all(durbanids_in_injtn_vl %in% durban0$RunID))


##  Save -------------------------------------------------------------------

#  confirm each `sampleid` is unique
stopifnot(anyDuplicated(injtn$sampleid) == 0)

save(injtn, file= fn$o$injn)
