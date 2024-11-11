This has some information and description about the code here, which generated
primary results of the project below. As this document is written in Markdown
(`.md`)., it is displayed better with a Markdown viewer.

## Project info

* Title : Quantifying the dynamics of the blood plasma proteome during acute HIV-1 infection
* NBIS Support #4964 -> #5800
* Principal Investigator: Joakim Esbjörnsson (<joakim.esbjornsson@med.lu.se>)
* Researcher : Jamirah Nazziwa (<jamirah.nazziwa@med.lu.se>)
* Organization : LU 
* Subject : Proteomics
* NBIS experts : Eva Freyhult (<eva.freyhult@nbis.se>) and Mun-Gwan Hong (<mungwan.hong@nbis.se>)

--------------------------------------------------------------------------------

## General info

### File name

The name of each script file often reflects what it does. It follows this format.
`<what it does>-<sample set ID>-<details>.[version].<extension>`


### Script header

Every script has a **header** at the top of it. It has a description about the
script as well as basic info. All **input** and **output** file names are listed
in the header.

### Neighbor folders

Those input and output files are supposed to be stored in the folders described
below. Every script is written to be executed on this current folder (`code`)
accessing those folders using relative paths.

```
project/
 ├── code/              <The current folder>
 │
 ├── data/               Data 
 │   ├── raw_external/     from external sources, e.g. public databases
 │   ├── raw_internal/     raw input data, obtained internally.
 │   └── interim/          intermediate data
 ├── reports/            Reports
 └── results/            Main output
     ├── figures/
     └── tables/
```

### R and Rmd scripts

All R **functions** except those from 
`R-core` (e.g. `base`, `stats`, `utils`, ...), 
`tidyverse` (e.g. `dplyr`, `ggplot2`, `purrr`, ...),
`kabelExtra` and `janitor` packages
are primarily called with double colons (`::`) to specify source package
clearly.  The R code followed [the `tidyverse` **R style
guide**](https://style.tidyverse.org/syntax.html).

--------------------------------------------------------------------------------

## Restore environment

Different software environment, e.g. different versions of R packages, from the
one during the development of the code, can yield unexpected error message or
generate dissimilar analysis results. To make the results reproducible, the
software environment around the code are saved. It can be restored ahead of
using the code following the steps below. Here, the process was facilitated by
`Conda` management system. The environment handling software `conda` is assumed
to be pre-installed. Please refer to [Conda](https://docs.conda.io) for
questions about the installation.

1. Open Terminal (Mac / Linux) or Command line (Windows)

2. Navigate to the directory where this README.md file is located using the `cd`
command.

3. Run the commands below. It creates the same Conda environment and activate it

    ```sh
    conda env create -n nbis5800 -f environment.yml
    conda activate nbis5800
    Rscript -e "remotes::install_github('mikabr/ggpirate@master')"
    ```

The Conda environment was stored in this file.

* `environment.yml`


--------------------------------------------------------------------------------

## Files

### Master script files

Run these lines below in an R console, which will generate all intermediate data
files for data analysis and
the QC report. Please make sure the input files are ready in expected
folders. The required files are listed below the code lines.  

```R
source("master_data_generation.R")   # data generation
```
    ../data/raw_internal/2020-10-23/Patient_Matched_sampleIDs_updated_AH_02042020 (003).xlsx
                         2020-10-23/05. Sample_infomation_Durban.xlsx
                         2021-01-11/symptoms_to_Mun-Gwan_11012021.csv
                         2022-02-01/depleted_20220201.txt.gz
                                    neat_20220201.txt.gz 
                         2021-11-19/NeatPlasma_MS_run_conditions.xlsx

The dates in the folder names above indicate when the files were transferred to
NBIS. If multiple versions of one file were transferred, please use the one
delivered on that day. 

### Data generation

#### `master_data_generation.R`

The master file for R data generation. This executes following scripts in the
proper order. 

* `gen-s1-viral_load.v01.R` : Read viral load data
* `gen-s1-clinc.v01.R` : Read clinical info data
* `gen-s1_neat_Spectronaut.v01.R` : Read and parse the output of Spectronaut. Attach tables about injection.
* `gen-s1_depl_Spectronaut.v01.R` : Read depleted plasma data
* `gen-s1_Spectronaut.v02.R` : Preprocessing of the Spectronaut data


### Report writing

* `report-s1_neat_Spectronaut.01-QC.Rmd` : About QC and preprocessing of Spectronaut data

