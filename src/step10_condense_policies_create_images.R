# Process policy set data and create barplots representing decision variables
# Script modified from version provided by Reclamation in June 2022
# Modified by E Stark August 2023, then M Pernat August 2023

rm(list=ls())
options( java.parameters = "-Xmx4g" )
options(scipen = 999)

source('src/utils/library.R')

src_dir <- 'borg_directories/'
borg_experiments <- list.dirs(
  path = src_dir, 
  full.names = FALSE, 
  recursive = FALSE
)

archive_file_name <- '/borg_run/Archive.txt'

n_powell_tiers <- 5
n_mead_tiers <- 8

for (experiment in borg_experiments){
  
  # Create a folder for the experiment in 'output'
  output_dir <- paste0('output/borg_experiments_analysis/', experiment)
  dir.create(output_dir, recursive = TRUE)
  
  # Path of the actual borg experiment directory
  borg_directory <- paste0(src_dir, experiment)
  
  # Get experiment's xml file
  files <- list.files(borg_directory)
  xml_file_name <- files[grep("\\.xml$", files)]
  
  # Condense policies for plotting
  condensed_archive <- condense_policies(
    borg_directory = borg_directory, 
    archive_file_name = archive_file_name,
    xml_file_name = xml_file_name,
    n_powell_tiers = n_powell_tiers,
    n_mead_tiers = n_mead_tiers
  )
  
  # Save the condensed archive
  write.csv(
    x = condensed_archive, 
    file = paste0(output_dir, '/condensed_archive.csv'),
    row.names = FALSE
  )
  
  # Get DV indices of condensed archive
  dv_indices <- get_dv_indices(
    condensed = TRUE,
    col_names = colnames(condensed_archive),
    n_powell_tiers = n_powell_tiers,
    n_mead_tiers = n_mead_tiers
  )
  
  # Create policy images. 
  create_policy_images(
    output_dir = output_dir,
    archive_df = condensed_archive,
    dv_indices
  )
  
  # Create Mead policy heatmap
  mead_heatmap_matrix <- MeadHeatmapMatrix(
    condensed_archive = condensed_archive,
    dv_indices = dv_indices,
    max_tiers = n_mead_tiers,
    max_elev = 1200,
    min_elev = 895,
    disc_length = 5
  )
  
  # Optional, order policies (columns) by starting shortage elevation
  mead_heatmap_matrix <- SortByStartingShortageElevation(
    mead_heatmap_matrix = mead_heatmap_matrix,
    max_shortage_elev = 1145
  )
  
  # Draw the heatmap
  ## Only keep half of the elevation values for tick labels to reduce clutter
  rownames(mead_heatmap_matrix) <- ifelse(
    seq_along(rownames(mead_heatmap_matrix)) %% 2 == 0,
    "",
    rownames(mead_heatmap_matrix)
  )
  
  heatmap <- Heatmap(
    matrix = mead_heatmap_matrix,
    col = GetMeadColorScheme(),
    cluster_rows = FALSE, 
    show_column_names = FALSE, 
    show_column_dend = FALSE,
    column_title = experiment, 
    name = "Shortage Volume \n(kaf)", 
    border = "black", 
    cluster_columns = FALSE,
    row_names_side = "left"
  )
  
  png(paste0(output_dir, '/mead_policy_heatmap.png'), width=800, height=600)
  draw(heatmap, heatmap_legend_side = "right")
  dev.off()
}




