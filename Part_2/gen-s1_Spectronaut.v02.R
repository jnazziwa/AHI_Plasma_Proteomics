# -----------------------------------------------------------------------------#
# Quality control of Spectronaut data of both neat and depleted samples
# This is created primarily from `report-s1_Spectronaut.01-QC.Rmd` by `knitr::purl`
# -----------------------------------------------------------------------------#
# initiated on 2021-08-04
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())
stopifnot(getRversion() >= "4.1.0")   # R base pipe
# Please note some functions are called directly using '::' or ':::'
library(tidyverse)

#----- File names --------------------------------------------------------------

fn <- list(
  i = list(   # input
    nt01 = "../data/s1_neat_Spectronaut.v01.RData",
    dp01 = "../data/s1_depl_Spectronaut.v01.RData"
  ),
  o = list(   # output
    interim = "../data/interim/Spectronaut_after_trimming.RData",
    sp02 = "../data/s1_Spectronaut.v02.RData"
  )
)
# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))

#----- MAIN --------------------------------------------------------------------

## Load data
# Load and combine proteomic data
e <- lapply(
  setNames(c("nt01", "dp01"), c("Neat", "Depleted")), 
  \(.x) {
    load(fn$i[[.x]], e <- new.env())
    e
  }
)
proteins <- imap_dfr(e, ~mutate(.x$proteins, prep_type = .y)) |> 
  mutate(protein_type_id = paste(id, prep_type, sep = "-")) 
spct <- imap_dfr(e, ~mutate(.x$spct, prep_type = .y)) |> 
  mutate(protein_type_id = paste(protein_id, prep_type, sep = "-")) |> 
  separate(sampleid, c("subjid", "visit_nr"), sep = "-", remove = FALSE)
prep_types <- set_names(unique(spct$prep_type))
rm(e)

#' Transform a long data frame to a matrix. Row names were filled with given
#' values in `row`.
#'
#' @param x a longitudinal data frame that has those three columns given by
#' `row`, `col` and `value`
#' @param row,col which column will be the rows/columns of the output
#' @param value the column of the values of the output matrix
#'
#' @return a matrix
long_to_matrix <- function(x, 
                           row = sampleid, 
                           col = protein_type_id, 
                           value = LogIntensities) {
  pivot_wider(
    data = x,
    id_cols = {{row}}, 
    names_from = {{col}},
    values_from = {{value}}
  ) |> 
    column_to_rownames(as.character(enexpr(row))) |> 
    as.matrix()
}

#' Transform to a list of matrices, which were divided by `prep_type`
#'
#' @param x a data frame
#' @return a list of matrices
to_matrix_by_prep_type <- function(x) {
  split(x, x$prep_type) |>     # by prep_type
    lapply(long_to_matrix)
}

## Remove the orphan sample
spct <- spct |> filter(sampleid != "194532-T1")

## * Replace below zero *
spct <- spct |> 
  mutate(LogIntensities = if_else(LogIntensities < 0, 0, LogIntensities))

#' Compute the proportion of a condition, together with other summary statistics
#'
#' @param mat a input matrix
#' @param cond_fun a function that returns a logical vector per row/column. The
#'   proportion of TRUE is stored in the `prop` of the output.
#' @param aggre_funs a list of functions with names, each of which compute a
#'   summary statistic per row/column, e.g. median. Every function should accept
#'   the argument `na.rm`. As the default, `na.rm` is set to TRUE.
#' @param wises a vector of two characters, which tell what rows/columns
#'   contain. The `wise` column of the output will get them indicating
#'   directions.
#' @return a tibble that has `wise`, `id`, `prop` and the names of the
#'   functions. The `prop` has the proportions of the values satisfied the
#'   condition given by `cond_fun`. Each column with a function name contains 
#'   the output of the function.
#' @importFrom rlang as_function
cond_prop <- function(mat, cond_fun, aggre_funs = NULL, wises = c("Sample", "Protein")) {
  # names are used for column names
  nms <- names(aggre_funs)
  if (!is.null(aggre_funs) && is.null(nms) || any(nms == "")) {
    stop("The names of `aggre_funs` must be given.")
  }
  
  lapply(
    1L:2L,    # row-wise , column-wise
    function(ii) {
      tibble(
        wise = case_when(
          ii == 1L ~ wises[1], 
          ii == 2L ~ wises[2],
          TRUE ~ ""
        ),
        id = dimnames(mat)[[ii]],
        prop = apply(mat, ii, \(x) {
          rlang::as_function(cond_fun)(x) |>  # purrr lambda function
            mean(na.rm = TRUE)
        }),
      ) |> 
        add_column(
          map_dfc(aggre_funs, ~apply(mat, ii, rlang::as_function(.x)))
        )
    }
  ) |> 
    bind_rows() |> 
    mutate(wise = factor(wise, levels = wises))   # fix the order for plots
}

cond_prop_by_prep_type <- function(x, cond_fun, aggre_funs = NULL) {
  to_matrix_by_prep_type(x) |> 
    # compute missing proportion and other statistics
    map(~ cond_prop(.x, cond_fun, aggre_funs)) |>
    bind_rows(.id = "prep_type")   # longitudinal table
}

## Missing proportion
sumstat <- list(
  Median = ~median(.x, na.rm = TRUE) 
  # MAD = ~mad(.x, na.rm = TRUE)
) 

# contains missing proportion together with median, MAD
mssngs <- cond_prop_by_prep_type(spct, is.na, sumstat) |> 
  expss::apply_labels(prop = "The proprotion of missing values")

## * Remove SAMPLES with too many missing values *
too_many_mssng_smpl <- mssngs |> 
  filter(prop > 0.7 & wise == "Sample") |> 
  arrange(desc(prop))
spct <- spct |> 
  anti_join(too_many_mssng_smpl, by = c("sampleid" = "id", "prep_type"))

## Missing proportions updated
mssngs_upd <- cond_prop_by_prep_type(spct, is.na, sumstat) |> 
  expss::apply_labels(prop = "The proprotion of missing values")

## * Remove PROTEINS with too many missing values *
too_many_mssng_protein <- mssngs_upd |> 
  filter(prop > 0.8 & wise == "Protein") |> 
  arrange(desc(prop))
spct <- spct |> 
  anti_join(too_many_mssng_protein, by = c(protein_type_id = "id"))

## Intermediate status save
save(spct, file = fn$o$interim)

## Imputation
### A modified version of Christofer Karlsson's approach
set.seed(20220202)

spct_imp <- lapply(
  split(spct, spct$prep_type),
  function(one_type) {
    # Non-log-transformed intensities
    spct_not_log <- one_type |> 
      pivot_wider(sampleid, protein_type_id, values_from = LogIntensities) |>
      pivot_longer(-sampleid, "protein_type_id", values_drop_na = FALSE) |>
      mutate(value = 2^value)
    
    # determine max impute value per protein
    max.impute.value <- spct_not_log |> 
      group_by(protein_type_id) |>
      summarise(max_imp = min(value, na.rm = TRUE), .groups = "drop")
    
    # join the max impute values to the prot.df
    # impute missing values with runif
    prot.df.imputed <- spct_not_log |>
      left_join(max.impute.value, by = "protein_type_id") |>
      mutate(
        Imputed = if_else(is.na(value), TRUE, FALSE),
        r_unif = sapply(max_imp, function(x) {
          runif(1, min = 1, max = x)
        }),
        value = if_else(is.na(value), r_unif, value)
      )
    
    # Store imputed data
    prot.df.imputed |>
      mutate(LogIntensities = log2(value)) |>
      select(sampleid, protein_type_id, LogIntensities)
  }
) |> 
  bind_rows(.id = "prep_type")

## LOESS normalization
# Normalized by Cyclic Loess
spct_proc <- spct_imp |>
  # per preparation type
  to_matrix_by_prep_type() |>
  map_dfr(
    ~t(.x) |>
      limma::normalizeCyclicLoess() |>
      as_tibble(rownames = "protein_type_id") |>
      pivot_longer(-protein_type_id, "sampleid", values_to = "LogIntensities")
  )


## Save
proteins <- proteins %>%
  filter(protein_type_id %in% spct_proc$protein_type_id)

save(spct_proc, proteins, file= fn$o$sp02)

