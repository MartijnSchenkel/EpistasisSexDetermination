# ================== #
# ---- METADATA ---- #
# ================== #
Data repository associated with:

Schenkel MA, Beukeboom LW & Pen I (2020). Epistatic interactions between sex chromosomes and 
autosomes can affect the stability of sex determination systems. [JOURNAL, VOLUME, PAGES, DOI]

# ======== #
# CONCENTS #
# ======== #

The repository contains the following folders, each of which contains a file readme.txt 
describing the files within:

/data_analysis_visualisation/

  Scripts for data collection, analysis, and visualisation

/figures

  Source figure files for Figures 1-3 and Supplementary Figures 1-2.

/generate_run_job_scripts/

  Scripts to generate and run jobscripts for simulations. Requires generating source files using
  files in /generate_source_scripts/

/generate_source_scripts/

  Scripts to generate the source code for all the different simulations for the manuscript. Requires
  seed generators to be compiled (found in /seed_generators/). Note that these are Bash scripts that
  generate .R files containing the actual source code that is run in R. Ergo the source code is not
  actually in an .R file itself, but rather is generated.

/secondary_data/
  
  Contains the secondary data file used to perform the data analyses. Included for the sake of ease 
  in case of future re-analysis.

/seed_generators/

  Source code for the seed generators used.


# ========================= #
# NOTE ON GENE NOMENCLATURE # 
# ========================= #

The files used here assume a slightly different nomenclature with regard to gene names as referred
to in the manuscript. In short, gene names are as follows (manuscript name = script name):

* Y   = YM
* A   = AM
* W   = FD
* SAA = SAI
* SAW = SAII

