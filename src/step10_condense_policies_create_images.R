# Process policy set data and create barplots representing decision variables
# Script modified from version provided by Reclamation in June 2022
# Modified by E Stark August 2023, then M Pernat August 2023

rm(list=ls())
options( java.parameters = "-Xmx4g" )
options(scipen = 999)

require(devtools)
library(tidyr)
library(patchwork)
source('src/utils/library.R')

borg_experiment_directories <- list.dirs(
  path = 'borg_directories/', 
  full.names = TRUE, 
  recursive = FALSE
)
archive_file_name <- '/borg_run/Archive.txt'

for (directory in borg_experiment_directories){
  
  # Get policies into condensed form for plotting
  condensed_archive <- condense_policies(
    directory = directory, 
    archive_file_name = archive_file_name,
    n_powell_tiers = 5,
    n_mead_tiers = 8
  )
  
  # Save the archived as condensed policies
  write.csv(
    x = condensed_archive, 
    file = paste0(directory, '//condensed_archive.csv'),
    row.names = FALSE
  )
  
  # Create policy images. Saved in the borg directory.
  create_policy_images(
    directory = directory,
    archive_df = condensed_archive
  )
  
}




