# -----------------------------------------------------------------------------#
# This can be run in R to get all intermediate data files as below.
#    source("master_data_generation.R")
# It runs all other R scripts that generate data files for data analysis. 
# -----------------------------------------------------------------------------#
# initiated on 2021-01-25
# authors : Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

library(dplyr)

#  Read the table and create RData file of viral load
source("gen-s1-viral_load.v01.R")

#  Read the table and create RData file of clinical info
source("gen-s1-clinc.v01.R")

# * "Neat" samples without any depletion
# 1. Read and parse the output of Spectronaut
# 2. Read and parse the tables about injection
# 3. make `durban` clean according to info (meeting note : 2020-10-26)
# 4. make `injtn` tidy
source("gen-s1_neat_Spectronaut.v01.R")

# * Proteomic data of depleted plasma samples
# 1. Read and parse the output of Spectronaut
# 2. Make a list of UniProt ID for conversion table
source("gen-s1_depl_Spectronaut.v01.R")


# Quality control of Spectraonaut proteomics data of both sample types
#   1) Too many missing
#   2) Imputation and Normalization
#   *  Most parts were extracted from `report-s1_Spectronaut.01-QC.Rmd`
source("gen-s1_Spectronaut.v02.R")


