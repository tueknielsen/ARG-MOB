# ARG-MOB
This repository provides documentation for the analysis in the preprint manuscript "Mobilization of antibiotic resistance genes differ by resistance mechanism" by Tue Kjærgaard Nielsen, Patrick Denis Browne, and Lars Hestbjerg Hansen. https://www.biorxiv.org/content/10.1101/2021.01.10.426126v1

"ARG_MOB_v0_5.sh" is the main script and contains the most information. Please edit paths and variables in the top of the script to match your paths before running the script. 

It is required that you also download the file "ARGMOB_dbs.tar.gz", before running the script. This file contains all necessary data files to run the analyses. For now, this file is only provided upon request to tkn@plen.ku.dk. 
The database files in ARGMOB_dbs.tar.gz include:
  - RefSeq complete bacterial genomes downloaded on Dec. 12 2019.
  - CARD database (v3.0.7). This database has not been modified from the original version. The original source of the database is McMaster University at https://card.mcmaster.ca/. The CARD database made available in ARGMOB_dbs.tar.gz is not in affiliation with or endorsed by McMaster University.
  - IS transposase protein sequences from ISfinder as implemented in Prokka (v1.14.6) under the GPLv3 license.

For the analysis, databases may be provided in a reformatted version. Please see ARG_MOB_v0_5.sh for notes on how edits were made. 

The "ARG_MOB_v0_5.sh" script will take a long time to finish. When it is done, the file "argmob_data.tar.gz" will be generated. This file is required input for the Rmarkdown script "ARGMOB_clean.rmd" which will perform statistics and generate plots. Before running the Rmarkdown script, please make sure you have the required packages installed (listed in the Rmarkdown document) and change the working directory to a suitable location. It is recommended to Knit the Rmarkdown document to a html output file, as this will allow for interactive tables and figures. 

