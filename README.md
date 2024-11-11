Welcome to the README for the GitHub repository housing the code used to generate primary results for the project detailed below. 
This document serves as a guide on how to replicate these results, providing insights into the code structure and necessary steps.
It's formatted using Markdown (.md), which renders best with a Markdown viewer.

**Project Title**: Quantifying the dynamics of the blood plasma proteome during acute HIV-1 infection.


--------------------------------------------------------------------------------

## Meet the Bioinformatics Team
This project was supported by NBIS and BIOMS and was executed by a team of four skilled Bioinformaticians/Statisticians: 

**Dr. Mun-Gwan Hong, NBIS**: Responsible for Quality Control and normalization of Spectronaut data in the project (*NBIS PROJECT ID: 4964).
**Dr. Eva Freyhult, NBIS**: Responsible for statistical and downstream analysis of blood plasma proteome dynamics during acute HIV-1 infection.(*NBIS PROJECT ID: 5800)
**Filip Årman, BIOMS**: Managed MS data upload and reviewed all QC steps.
**Dr. Jamirah Nazziwa, LU**: Responsible for data visualization derived from the project.

--------------------------------------------------------------------------------

## Project Details

* Title : Quantifying the dynamics of the blood plasma proteome during acute HIV-1 infection
* NBIS Support #4964 -> #5800
* Principal Investigator: Joakim Esbjörnsson (<joakim.esbjornsson@med.lu.se>)
* Researcher : Jamirah Nazziwa (<jamirah.nazziwa@med.lu.se>)
* Organization : LU 
* Subject : Proteomics
* NBIS experts : Eva Freyhult (<eva.freyhult@nbis.se>) and Mun-Gwan Hong (<mungwan.hong@nbis.se>)

--------------------------------------------------------------------------------

## General information

### File Naming Convention

R Scripts are named according to their functionality, sample set ID, and additional details.
It follows this format.
`<what it does>-<sample set ID>-<details>.[version].<extension>`


### Script structure

Each script begins with a **header** containing a description and basic info, including **input** and **output** file names

### Folder Structure

The scripts are designed to work within a specific folder structure, facilitating relative path access to input and output files.
Every script is written to be executed on this current folder (`code`)accessing those folders using relative paths.



```
project/            #Choose a name for your project
 ├── code/              #<The current folder: the directory that includes this README.md>
 │   ├── Renv           #contains all packages and versions │
 ├── data/               #Data 
 │   ├── raw_external/     #from external sources, e.g. public databases
 │   ├── raw_internal/     #raw input data, obtained internally.
 │   ├── intermediate/     #intermediate data for Data analysis.
 │   └── interim/          #intermediate data for QC analysis
 ├── reports_aux/            #Reports
 └── results/            #Main output
     ├── figures/
     └── tables/
```

If these folders aren't set up, instructions can be found in the Appendix.


--------------------------------------------------------------------------------
# START HERE

The analysis pipeline was divided into 4 parts .i.e.

Part 1: Estimating relative protein abundances from peptide quantification
Part 2: QC assessment Project #SMS-4964-19-hiv
Part 3: Data analysis with Project #SMS-5800-HIV
Part 4: Visualisation for publication

--------------------------------------------------------------------------------
# PART 1: Estimating relative protein abundances from peptide quantification
1. Download the spectronaunt files uploaded at ProteomeXchange via the PRIDE database 

Depleted plasma file: "20220201_073608_DEPplasmaFAIMS-HIV-Durban+IAVI_CombinedLibrarySearch_210608_Report_iq.xls"
Neat Plasma file: "20220201_110528_ACUTEHIV_Report_iq-xls.xls"

ProteomeXchange with identifier PXD042850

2. Transfer these xls files (about 6GB each file) to the Part_1 folder
3. Make Part_1 folder your working directory in R Studio
4. Open R script maxlfq.r (inside Part_1 folder) in R-studio and run it.
5. The script above will produce estimated relative protein abundances for all proteins.
The resulting files from this part were  "neat_20220201.txt" and "depleted_20220201.txt" 
6. These files are the input files in Part 2 on QC assessment

--------------------------------------------------------------------------------
# PART 2: Reproducing Results for QC assessment Project #SMS-4964-19-hiv

These scripts were developed in R Version 4.1 and version 4.2 for visuliasation

To reproduce the results, follow these steps:

## Step 1: Create an Rstudio project
Set the working directory to the `code/` folder

## Step 2: Prepare Directories
Ensure input and output directories match the specified structure

## Step 3: Copy code (files with .R extension) from SMS-4964 (Inside Analysis >Git_repo> QC> SMS-4964-19-hiv)

Please Note the files have alrady been transfered in the code folder, but if you want to start from scratch, then move these files yourself

Transfer all `.R` files from the provided source to`<project>/code/`.
Transfer Renv folder to `<project>/code/`.
Transfer Renv lock file to `<project>/code/`. 

### R and Rmd scripts

All R **functions** except those from 
`R-core` (e.g. `base`, `stats`, `utils`, ...), 
`tidyverse` (e.g. `dplyr`, `ggplot2`, `purrr`, ...),
`kabelExtra` and `janitor` packages
are primarily called with double colons (`::`) to specify source package clearly.  
The R code followed [the `tidyverse` **R style guide**](https://style.tidyverse.org/syntax.html).


## Step 4: Restore environment

Different software environment, e.g. different versions of R packages, from the one during the development of the code, can yield unexpected error message or generate dissimilar analysis results. 
To make the results reproducible, the software environment around the code are saved. It can be restored ahead of using the code following the steps below. 
Here, the process was facilitated by `Conda` management system. The environment handling software `conda` is assumedto be pre-installed. Please refer to [Conda](https://docs.conda.io) forquestions about the installation.

### Conda 

1. Open Terminal (Mac / Linux) or Command line (Windows)

2. Navigate to the directory where this README.md file is located using the `cd`
command.

3. Run the commands below. It creates the same Conda environment and activate it

If you prefer, [Mamba](https://mamba.readthedocs.io/) can be used instead of Conda.
    ```sh
    conda env create -n nbis5800 -f environment.yml
    conda activate nbis5800
    Rscript -e "remotes::install_github('mikabr/ggpirate@master')"
    ```

The Conda environment was stored in this file.

* `environment.yml`

### R (renv)

Execute these commands below in the `R`. 
These will restore the same R environment using `renv` R package and reproduce the primary results including the final report.

```R
renv::settings$external.libraries(.libPaths())
renv::activate()
renv::restore()
```
More information about the R environment handling package,`renv` can be found at  ["introduction to `renv`"](https://rstudio.github.io/renv/articles/renv.html).

## Step 5: Execute Code for Analysis

Run the master script files in R console to generate intermediate data files and QC reports.  

```R
source("master_data_generation.R")   # data generation
```


--------------------------------------------------------------------------------

## REQUIRED Files

Ensure the availability of specific files in the designated folders for successful execution.

    ../data/raw_internal/2020-10-23/Patient_Matched_sampleIDs_updated_AH_02042020 (003).xlsx
                         2020-10-23/05. Sample_infomation_Durban.xlsx
                         2021-01-11/symptoms_to_Mun-Gwan_11012021.csv
                         2021-04-03/Protocol_C_with_ART_start_date_for_JN_23042021.xlsx
                         2021-04-23/Protocol_C_with_ART_start_date_for_JN_23042021.xlsx
                         2021-06-28/protein_dep_20210609.txt.gz
                         2021-11-25/Copy of Depleted_plasma_samples_DIAresultSummary_plus depletionbatch_19Nov21.xlsx
                         2021-11-25/NeatPlasma_MS_run_conditions.xlsx
                         2022-02-01/depleted_20220201.txt.gz
                         2022-02-01/neat_20220201.txt.gz 
                         2021-11-19/NeatPlasma_MS_run_conditions.xlsx
						 CD4_CD8_Proteomics_08062021.xlsx
						 DEP_processing_date.xlsx
						 Durban_ART_startdate.xlsx
						 Durban_CD4_VL_HLA_KIR.xlsx
						 IAVI_CD4_VL_HLA_KIR.xlsx

The dates in the folder names above indicate when the files were transferred to NBIS. 
If multiple versions of one file were transferred, please use the one delivered on that day. 

### Data generation

#### `master_data_generation.R`

The master file for R data generation. This executes following scripts in the
proper order. 

* `gen-s1-viral_load.v01.R` : Read viral load data
* `gen-s1-clinc.v01.R` : Read clinical info data
* `gen-s1_neat_Spectronaut.v01.R` : Read and parse the output of Spectronaut. Attach tables about injection.
* `gen-s1_depl_Spectronaut.v01.R` : Read depleted plasma data
* `gen-s1_Spectronaut.v02.R` : Preprocessing of the Spectronaut data


### Step 6: Report writing for QC reporting
All .Rmd should be in the REPORT_AUX folder

* `report-s1_neat_Spectronaut.01-QC.Rmd` : About QC and preprocessing of Spectronaut data


```R
# QC report writing
rmarkdown::render("report_aux/report-s1_neat_Spectronaut.01-QC.Rmd", output_dir = "../report_aux/")     
```

Output file: report-s1_Spectronaut.01-QC.html Size: 7,858KB
* This takes about 30 mins on a 32GB RAM computer
--------------------------------------------------------------------------------

# PART 3: Reproducing Results for Data analysis with PROJECT #SMS-5800-HIV
This part depends on outputs from QC assessment Project **#SMS-4964-19-hiv**. 
These files were created after **successful execution** up to Part_2 above and they were written in the Data folder  

*s1-viral_load.v01.RData*
*s1-clinic.v01.RData*
*s1_neat_injection.v01.RData
*s1_neat_Spectronaut.v01.RData
*s1_depl_Spectronaut.v01.RData
*s1_Spectronaut.v02.RData


## Step 1: Move all.RMD files in project SMS-5800-HIV in the report_aux FOLDER

data_inj_v0.Rmd
data_neatdeploverlap_v0.Rmd
data_readneatdepl_v0.Rmd
data_readneatdepl_v1.Rmd
data_readneatdepl_v2.Rmd
data_readneatdepl_v3.Rmd
data_readneatdepl_v4.Rmd
evalclusterCD4abs.Rmd
evalclusterVL112m.Rmd
evalclusterXXXid.Rmd
header.tex
NBIS5800neatdepl.Rmd
splineclusterCD4abs.Rmd
splineclusterVL.Rmd
splineclusterXXXid.Rmd

The header.tex file is used to create logos for the report. Edit the file directory of these logos inlines 22 and 23 of this tex file or delete them.

## Step 2: Perform data analysis steps
This step produces a report of about 107 pages that includes all the analysis done to quantiate the HIV-1 plasma proteome. Run these line below in an R console, which will generate the report.

```R
# Data analysis report writing
rmarkdown::render("report_aux/NBIS5800neatdepl.Rmd", output_dir = "../report_aux/")     
```


# FOR THE REVIEW PROCESS: WE UPDATED THE ANALYSIS BASED ON REVIEWER COMMENTS
Use the qmd file instead---NBIS5800review.qmd.
```R
# Data analysis report writing
install.packages("quarto")
quarto::quarto_render("report_aux/NBIS5800review.qmd")     
```


* This takes about 8 hours on a 32GB RAM computer


### Output File
NBIS5800neatdepl.pdf : size 24,983KB

--------------------------------------------------------------------------------

# PART 4: Reproducing Visualisation used in Manuscript

All files and results needed for the plots are uploaded visualisation folder
The generated plots are saved as pdf, these can be excluded when reproducing the plots
The rscripts are labeled depending on the outcome (ARS, VL, dynanamics) as outlined in the manuscript
Input files were generated from Part 3.

Restore the R environment as described in part 1, using Part_3
R scripts are labelled based on the different sections in results


--------------------------------------------------------------------------------

# Appendix
Below are some useful commands to create directories and set up the project structure. 
Follow these commands in your terminal or command line interface.

That's it! You're all set to reproduce the analysis results. If you encounter any issues or have questions, feel free to reach out to the project team for assistance.
## Useful commands

Here are some potentially useful commands.
Note that these are supposed to run on Terminal (Mac OS/Linux) or Command line (Windows), depending on the (operating) system you have.

### Directories

If those directories are not available yet, they can be created by the commands below.
Find an approriate place for the project.
Move there first, using `cd` command.

#### Mac/Linux

    mkdir nbis_support &&
      cd nbis_support &&
      mkdir -p code/ \
               data/{external,raw_internal,processed,intermediate}/ \
               reports/ \
               results/{figures,tables}/

#### Windows

    mkdir nbis_support &&^
      cd nbis_support &&^
      md code\ data\external\ data\raw_internal\^
         data\processed\ data\intermediate\^
         reports\ results\figures\ results\tables\

Then, copy all the files and subdirectories in here to the newly created directory `code` under `nbis_support`. 
