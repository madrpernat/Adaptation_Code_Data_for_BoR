require(devtools)
library(ComplexHeatmap)
library(tidyr)
library(patchwork)
library(xml2)
library(plyr)
library(kohonen)
library(lhs)
library(stats)
library(aweSOM)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(data.table)
library(plotrix)
library(readxl)
library(circlize)

############################# SOM HELPER FUNCTIONS #############################
##################### Coded by Nathan Bonham - March 2022 ######################
################ (https://doi.org/10.1016/j.envsoft.2022.105491) ###############

# Function to transform continuous values to discrete categories
continuous_to_discrete <- function(unit_cube, probs, categories) {
  for (i in 1:(length(probs)-1)) {
    indices <- which(unit_cube > probs[i] & unit_cube <= probs[i+1])
    unit_cube[indices] <- categories[i]
  }
  return(unit_cube)
}

quantile2radius=function(fraction, x,y,shape){
  
  # I took this code from the source code of the supersom function in Kohonen 
  # package and turned it into a function to calculate the radius given a 
  # fraction of the max node-to-node distance
  
  grid = somgrid(x,y,shape)
  grid <- check.somgrid(grid)
  nhbrdist <- unit.distances(grid)
  radius = quantile(nhbrdist, fraction)
  
  return(radius)
  
}


check.somgrid=function (grd) # taken from Kohonen source code
{
  mywarn <- FALSE
  if (is.null(grd$toroidal)) {
    mywarn <- TRUE
    grd$toroidal <- FALSE
  }
  if (is.null(grd$neighbourhood.fct)) {
    mywarn <- TRUE
    grd$neighbourhood.fct <- factor("bubble", levels = c("bubble", 
                                                         "gaussian"))
  }
  if (mywarn) 
    warning("Added defaults for somgrid object - ", "you are probably using the somgrid function ", 
            "from the class library...")
  grd
}


unit.distances=function (grid, toroidal) # taken from Kohonen source code
{
  if (missing(toroidal)) 
    toroidal <- grid$toroidal
  if (!toroidal) {
    if (grid$topo == "hexagonal") {
      return(as.matrix(stats::dist(grid$pts)))
    }
    else {
      return(as.matrix(stats::dist(grid$pts, method = "maximum")))
    }
  }
  np <- nrow(grid$pts)
  maxdiffx <- grid$xdim/2
  maxdiffy <- max(grid$pts[, 2])/2
  result <- matrix(0, np, np)
  for (i in 1:(np - 1)) {
    for (j in (i + 1):np) {
      diffs <- abs(grid$pts[j, ] - grid$pts[i, ])
      if (diffs[1] > maxdiffx) 
        diffs[1] <- 2 * maxdiffx - diffs[1]
      if (diffs[2] > maxdiffy) 
        diffs[2] <- 2 * maxdiffy - diffs[2]
      if (grid$topo == "hexagonal") {
        result[i, j] <- sum(diffs^2)
      }
      else {
        result[i, j] <- max(diffs)
      }
    }
  }
  if (grid$topo == "hexagonal") {
    sqrt(result + t(result))
  }
  else {
    result + t(result)
  }
}


################################################################################

# This function parses the archive file assuming that columns are in the 
# following order: decision variables, objectives, constraints, metrics.
condense_policies <- function(
    
    borg_directory, 
    archive_file_name, 
    xml_file_name, 
    n_powell_tiers, 
    n_mead_tiers
    
  ){
  
  # Get objective, metric, and constraint info from xml file
  xml_file_path <- paste0(borg_directory, '/', xml_file_name)
  xml_info <- get_xml_info(xml_file_path)
  
  objectives <- xml_info$objectives
  metrics <- xml_info$metrics
  constraints <- xml_info$constraints
  
  n_objectives <- length(objectives)
  n_metrics <- length(metrics)
  n_constraints <- length(constraints)
  
  
  # Read in archive file
  archive_file_path = paste0(borg_directory, archive_file_name)
  archive_df <- read.table(archive_file_path, header=TRUE)
  

  
  # Get column indices for decision variables
  dv_indices <- get_dv_indices(
    condensed = FALSE,
    col_names = colnames(archive_df),
    n_powell_tiers = n_powell_tiers,
    n_mead_tiers = n_mead_tiers
  )
  ## Powell DVs
  PTierEl_idx <- dv_indices$PTierEl
  PRels_idx <- dv_indices$PRels
  MeadRefs_idx <- dv_indices$MeadRefs
  BalMaxOffset_idx <- dv_indices$BalMaxOffset
  BalMinOffset_idx <- dv_indices$BalMinOffset
  
  ## Mead DVs
  MeadSurplus_idx <- dv_indices$MeadSurplus
  MeadEl_idx <- dv_indices$MeadEl
  MeadV_idx <- dv_indices$MeadV
  
  # Get number of DV columns
  n_powell_columns <- length(
    c(PTierEl_idx, PRels_idx, MeadRefs_idx, BalMaxOffset_idx, BalMinOffset_idx)
  )
  n_mead_columns <- length(
    c(MeadSurplus_idx, MeadEl_idx, MeadV_idx)
  )
  n_dvs <- n_powell_columns + n_mead_columns
  
  # Drop the constraints columns
  if (n_constraints > 0){
    
    constraint_start <- n_dvs + n_objectives + 1
    constraint_end <- constraint_start + n_constraints - 1
    
    archive_df <- archive_df[, -c(constraint_start:constraint_end)]
  }

  
  # Create a Powell dataframe (PTiering) for each policy
  for (j in 1:nrow(archive_df)){
    
    PTierEl = t(archive_df[j, PTierEl_idx])
    PRels = t(archive_df[j, PRels_idx])
    MeadRefs = t(archive_df[j, MeadRefs_idx])
    BalMaxOffset = t(archive_df[j, BalMaxOffset_idx])
    BalMinOffset = t(archive_df[j, BalMinOffset_idx])
    
    PTiering = as.data.frame(cbind(
      PTierEl, 
      MeadRefs, 
      BalMaxOffset, 
      BalMinOffset, 
      PRels
    ))
    row.names(PTiering) = c(
      "PT1", 
      "PT2", 
      "PT3", 
      "PT4", 
      "PT5"
    )
    colnames(PTiering) = c(
      "PTierEl", 
      "MeadRefEl", 
      "BalMaxOffset", 
      "BalMinOffset", 
      "PRels"
    )
    
    # If there are multiple tiers at the same elevation, set their Mead 
    # Reference Elevations, Min Offset, Max Offset, and Primary Release values
    # (columns 2-5) to that of the row with the highest Primary Release volume.
    for (i in 1:(nrow(PTiering) - 1)){
    
      if (PTiering[i,1] == PTiering[i+1,1]){
        PTiering[i+1,2] = PTiering[i,2]
        PTiering[i+1,3] = PTiering[i,3]
        PTiering[i+1,4] = PTiering[i,4]
        PTiering[i+1,5] = PTiering[i,5]
      } else { 
        PTiering[i+1,2] = PTiering[i+1,2]
        PTiering[i+1,3] = PTiering[i+1,3]
        PTiering[i+1,4] = PTiering[i+1,4]
        PTiering[i+1,5] = PTiering[i+1,5]
      }
      
    }
    
    # Create bottom tier (if bottom row elev = 3370, lowest tier > 3370 becomes 
    # bottom tier; mead ref, primary rel w/ 99999; just keep elev and release 
    # range)
    if (PTiering[nrow(PTiering), 1] == 3370 && PTiering[1, 1] != 3370){
      
      next_min = which(
        PTiering[, 1] == min(PTiering$PTierEl[PTiering$PTierEl > 3370])
      )
      
      PTiering[next_min,"MeadRefEl"] = PTiering[next_min,"PRels"] = 99999999
      
    }
    
    # Make MeadRefEl = 99999999 if Min and Max offset are 0
    for (i in 1:(nrow(PTiering))){
      
      if (PTiering[i, 3] == 0 && PTiering[i, 4] == 0){
        
        PTiering[i, 2] = 99999999
        
      }
    }
    
    # If two adjacent tiers are the same except for their elevations, make them
    # one tier.
    for (i in 1:(nrow(PTiering) - 1)){
      
      for (k in (i + 1):nrow(PTiering)){
        
        if (all(PTiering[i, 2:5] == PTiering[k, 2:5])){
          
          PTiering[k, ] <- PTiering[i, ]
          
        } else{
          
          break
          
        }
      }
    }
    
    # If there are repeat rows, keep one of them - the others replace their pool 
    # elevation to 3370 and other values to 99999999. If any rows already have 
    # a pool elevation of 3370, replace their other values with 99999999 as well
    for (i in 1:(nrow(PTiering))){
      
      if (any(PTiering[i, 1] == PTiering[-i, 1]) | PTiering[i,1] == 3370){
        
        PTiering[i,1] = 3370
        PTiering[i,2] = 99999999
        PTiering[i,3] = 99999999
        PTiering[i,4] = 99999999
        PTiering[i,5] = 99999999
        
      } else { 
        
        PTiering[i,1] = PTiering[i,1]
        PTiering[i,2] = PTiering[i,2]
        PTiering[i,3] = PTiering[i,3]
        PTiering[i,4] = PTiering[i,4]
        PTiering[i,5] = PTiering[i,5]
        
      }
      
    }
    
    # Sort rows based on pool elevation (col 1)
    PTiering = PTiering[order(PTiering$PTierEl, decreasing = TRUE), ] 
    
    # 
    
    # Distribute columns of PTiering back into a single row and combine w/ 
    # un-condensed mead variables & objective values
    archive_df[j, PTierEl_idx] = t(PTiering$PTierEl)
    archive_df[j, MeadRefs_idx] = t(PTiering$MeadRefEl)
    archive_df[j, BalMaxOffset_idx] = t(PTiering$BalMaxOffset)
    archive_df[j, BalMinOffset_idx] = t(PTiering$BalMinOffset)
    archive_df[j, PRels_idx] = t(PTiering$PRels)
  }
  
  # MEAD
  
  ## Tier elevations
  short_elev = t(archive_df[, MeadEl_idx])
  
  ## Tier shortage volumes
  archive_df[, MeadV_idx] = apply(
    X = archive_df[, MeadV_idx], 
    MARGIN = 2, 
    FUN = function(x) as.numeric(x)
  )
  short_vol = archive_df[, MeadV_idx] / 1000
  short_vol = apply(short_vol, 1, rev) #reverse order of volumes
  
  # Initialize condensed dataframes
  compressed_vol = short_vol
  compressed_elev = short_elev
  
  # If repeating elevations, replace volume w/ highest. Replacing the elevations 
  # needs to come first b/c of how RW handles the variables
  for(i in 1:ncol(compressed_vol)){
    
    for(j in 1:nrow(compressed_vol)){
      
      if (any(compressed_elev[j,i] == compressed_elev[-j,i])){
        
        compressed_vol[j,i] = max(
          compressed_vol[which(compressed_elev[,i] == compressed_elev[j,i]), i]
        )
        
      }
    }
  }
  
  # If repeating volumes, replace elevation w/ highest
  
  for(i in 1:ncol(compressed_vol)){
    
    for(j in 1:nrow(compressed_vol)){
      
      if (any(compressed_vol[j,i] == compressed_vol[-j,i])){
        
        compressed_elev[j,i] = max(
          compressed_elev[which(compressed_vol[,i] == compressed_vol[j,i]), i]
        )
        
      }
    }
  }
  
  
  # Now replace all repeated rows w/ 895 & 99999999, then order descending and 
  # ascending to produce condensed tables
  
  for(i in 1:ncol(compressed_vol)){
    
    for(j in 2:nrow(compressed_vol)){
      
      if (any(compressed_vol[j,i] == compressed_vol[1:(j-1),i])){
        
        compressed_elev[j,i] = 895
        compressed_vol[j,i] = 99999999
      }
    }
  }
  
  # Replace any rows that have vol = 0 w/ 895 and 999999999 (to address 0s in 
  # first tiers)
  
  for(i in 1:ncol(compressed_vol)){
    
    for(j in 1:nrow(compressed_vol)){
      
      if (compressed_vol[j,i] == 0){
        
        compressed_elev[j,i] = 895
        compressed_vol[j,i] = 99999999
      }
    }
  }
  
  # Replace any rows that have elev = 895 w/ 999999999 (to address tiers that 
  # originally had volumes w/ elevation = 895)
  
  for(i in 1:ncol(compressed_vol)){
    
    for(j in 1:nrow(compressed_vol)){
      
      if (compressed_elev[j,i] == 895){
        
        compressed_vol[j,i] = 99999999
      }
    }
  }
  
  # Save row names, since they get wiped after 'apply' below
  row_names <- c(row.names(compressed_elev), rev(row.names(compressed_vol)))
  
  # Sort based on elevation only
  compressed_elev = apply(compressed_elev, 2, as.numeric)
  condensed_elev = apply(compressed_elev, 2, sort, decreasing=T)
  condensed_vol = apply(compressed_vol, 2, sort, decreasing=F)
  
  # Re-combine
  condensed_policies = as.data.frame(rbind(condensed_elev, condensed_vol))
  
  row.names(condensed_policies) = row_names
  
  colnames(condensed_policies) = c(1:ncol(condensed_policies))
  
  condensed_policies = t(condensed_policies)
  
  # Calculate elevations for surplus tier
  
  surplus_size = data.frame(archive_df[, MeadSurplus_idx])
  surplus_elev = 1200 - surplus_size
  
  # If surplus went away, replace 1200s w/ 99999999
  surplus_elev[surplus_elev == 1200] = 99999999
  colnames(surplus_elev) = "surplus_elev"
  
  # Add surplus elevation to condensed policy df
  condensed_policies_full = cbind(surplus_elev, condensed_policies)
  
  # NEW ARCHIVE_DF
  archive_df = cbind(
    condensed_policies_full,                    # Mead DVs
    archive_df[, PTierEl_idx],                  # Powell DVs
    archive_df[, PRels_idx],                    # Powell DVs
    archive_df[, MeadRefs_idx],                 # Powell DVs
    archive_df[, BalMaxOffset_idx],             # Powell DVs
    archive_df[, BalMinOffset_idx],             # Powell DVs
    archive_df[, (n_dvs + 1):ncol(archive_df)]  # Objectives and Metrics
  )
  
  col_names = sapply(
    X=colnames(archive_df), 
    FUN=column_mapping, 
    USE.NAMES = FALSE
  )
  
  names(archive_df) = col_names
  
  # Return archive with condensed policies
  archive_df
  
}


create_policy_images <- function(output_dir, condensed_archive, dv_indices){
  
  ## Mead DVs
  MeadSurplus_idx <- dv_indices$MeadSurplus
  MeadEl_idx <- dv_indices$MeadEl
  MeadV_idx <- dv_indices$MeadV
  
  #lastDV = which(colnames(condensed_archive) == 'T8V')
  
  #dvs = condensed_archive[, 1:lastDV]
  
  # Prepare for plotting. Add full pool and deadpool columns
  
  #elevation = dvs[1:9]
  #elevation = cbind(rep(1220, nrow(elevation)), elevation)
  #names(elevation) = c("Top", names(elevation[2:10]))
  #elevation$dead.pool=895
  
  elevation <- condensed_archive[, c(MeadSurplus_idx, MeadEl_idx)]
  elevation <- cbind(rep(1220, nrow(elevation)), elevation)
  names(elevation) <- c('Top', names(condensed_archive[, c(MeadSurplus_idx, MeadEl_idx)]))
  elevation$dead.pool <- 895
  
  # When surplus tier doesn't exist, need to replace value to make elev_delta work
  
  for(j in 1:nrow(elevation)){
    
    if (elevation[j, 2] == 99999999){
      
      elevation[j, 2] = 1220
    }
  }
  
  # Save max elevation for ggplotting later
  
  elev_max = unlist(apply(elevation[2:10], 1, max))
  
  # Initialize elev_delta
  
  elev_delta = elevation[, 1:10]
  
  # Difference between subsequent tier elevations, used for plotting
  
  for(i in 2:11){
    
    elev_delta[i - 1] = elevation[i - 1] - elevation[i]
    
  }
  
  volume = condensed_archive[, MeadV_idx]
  volume[volume == 99999999] = 0
  volume$dead.pool = 0
  vol_max = apply(volume, 1, max)
  
  # Identify actual tiers. Sum number of tiers. Want this for filtering in web app
  
  tier_bin = volume > 0
  
  nTiers = apply(tier_bin, 1, sum)
  
  # Create shortage volume labels, ex) V = 500 KAF
  
  volume_labs = matrix(NA, nrow = nrow(volume), ncol = ncol(volume))
  
  for (i in 1:ncol(volume_labs)){
    
    volume_labs[, i] = paste(volume[, i], ' KAF', sep='')
    
  }
  
  volume_labs[volume_labs == "0 KAF"] = NA
  volume_labs = data.frame(volume_labs)
  
  # Add 2 columns to volume df so it becomes the same number of columns as 
  # elev_delta (which has top and surplus tiers now)
  
  volume = cbind(rep(NA, nrow(volume)), rep(NA, nrow(volume)), volume)
  names(volume) = c("Top","Surplus", names(volume[3:11]))
  
  volume_labs = cbind(rep("Surplus", nrow(volume_labs)), 
                      rep("Normal", nrow(volume_labs)),
                      volume_labs)
  
  names(volume_labs) = names(volume)
  
  # Add policy ID to each data frame
  elevation$policy = as.character(1:nrow(condensed_archive))
  elev_delta$dead_pool = 895
  elev_delta$policy = as.character(1:nrow(condensed_archive))
  volume$policy = as.character(1:nrow(condensed_archive))
  volume_labs$policy = as.character(1:nrow(condensed_archive))
  
  piv_col = which(colnames(elev_delta)=='dead_pool')
  
  # Wide to long format for ggplot
  elevation = pivot_longer(elevation, cols=1:piv_col, names_to = 'Tier')
  elev_delta = pivot_longer(elev_delta, cols=1:piv_col, names_to = 'Tier') 
  
  # Need to reorder the factor levels such that dead pool is last
  
  elev_delta$Tier = factor(elev_delta$Tier, levels = c("Top",
                                                       "surplus_elev",
                                                       "T1e",
                                                       "T2e",
                                                       "T3e",
                                                       "T4e",
                                                       "T5e",
                                                       "T6e",
                                                       "T7e",
                                                       "T8e",
                                                       "dead_pool"))
  
  elev_delta = dplyr::rename(elev_delta, 'delta'='value')
  
  # Stacked bar plot order depends on the level order
  
  volume = pivot_longer(volume, cols=1:piv_col, names_to = 'Tier')
  volume$value = as.numeric(volume$value)
  volume_labs = pivot_longer(volume_labs, cols=1:piv_col, names_to = 'Tier')
  volume_labs = data.frame(volume_labs, elevation = elevation$value)
  
  # Remove label for surplus when it doesn't really exist (can start at row 2 b/c 
  # row 1 is Top)
  
  for (i in 2:nrow(volume_labs)){
    
    if(volume_labs[i, 2] == "Surplus" && volume_labs[i, 4] == 1220){
      
      volume_labs[i - 1, 3] = NA
    }
  }
  
  
  # Add color to volume df
  
  volume_col = as.data.frame(read.csv('data/vol_gradient.csv'))
  volume_col$color = as.character(volume_col$color)
  
  volume$color = ""
  volume$value = as.numeric(volume$value)
  
  for(i in 1:nrow(volume)){
    
    if(is.na(volume[i, 3]) || as.numeric(volume[i,3]) == 0){
      
      volume[i,4] = "#000000"
    }
    else if(as.numeric(volume[i, 3]) > 0){
      
      p = as.numeric(volume[i, 3])
      volume[i,4] = volume_col[which(volume_col$volume == p), 2]
      
    }
  }
  
  for(i in 1:nrow(volume)){
    
    if (volume[i, 2] == "Top"){
      
      volume[i, 4] = "#6d46a5"
    }
    
    if (volume[i, 2] == "Surplus"){
      
      volume[i, 4] = "#4659a5"
    }
    
    if (volume[i, 2] == "PartSurplus"){
      
      volume[i, 4] = "#bdbfce"
    }
    
  }
  
  df = data.frame(elev_delta,
                  v_lab = volume_labs$value,
                  v_col = volume$color,
                  elevation = volume_labs$elevation,
                  volume = volume$value)
  
  df$Tier = rep(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), 
                length(unique(df$policy)))
  
  df$TierName = rep(c('Surplus', 'Normal', '1', '2', '3', '4', '5', '6', '7', '8', 'Dead Pool'), 
                    length(unique(df$policy)))
  
  # Change to number and add zeros for image ordering in tableau
  
  df$policy = as.numeric(df$policy)
  df$policy = sprintf("%04d", df$policy)
  
  ##################################  POWELL  ##################################
  
  firstpDV = which(colnames(condensed_archive) == 'PT1e')
  lastpDV = which(colnames(condensed_archive) == 'MinOffset5')
  
  pdv = condensed_archive[, firstpDV:lastpDV]
  
  #rel_ranges.df=read.table("powell_release_range_table.txt", header = T, sep = )
  
  #################### Create data frame for stacked bar plots ###################
  
  # Translate dummy variables into NA
  
  pdv[pdv == 99999999] = NA
  
  # Prepare for plotting. Add equalization and deadpool columns.
  
  p.elevation = pdv[1:5]
  p.elevation = cbind(rep(3700, nrow(p.elevation)), p.elevation)
  names(p.elevation) = c("Equalization", names(p.elevation[2:ncol(p.elevation)]))
  p.elevation$dead.pool = 3370
  
  # Initialize elev_delta
  
  p.elev_delta = p.elevation[, 1:6]
  
  # Difference between tier elevations, used for plotting
  
  for(i in 2:ncol(p.elevation)){
    
    p.elev_delta[i - 1] = p.elevation[i - 1] - p.elevation[i]
  }
  
  # Mead reference elevations
  mead_ref = pdv[, 11:15]
  mead_ref = cbind(rep(NA, nrow(mead_ref)), mead_ref)
  names(mead_ref) = c("Equalization", names(mead_ref[2:ncol(mead_ref)]))
  
  # Start with min offset
  
  minoffsets = pdv[21:25] / 1000000
  minoffsets = cbind(rep(NA, nrow(minoffsets)), minoffsets)
  names(minoffsets) = c("Equalization", names(minoffsets[2:ncol(minoffsets)]))
  
  # Max offsets
  
  maxoffsets = pdv[16:20] / 1000000
  maxoffsets = cbind(rep(NA, nrow(maxoffsets)), maxoffsets)
  names(maxoffsets) = c("Equalization", names(maxoffsets[2:ncol(maxoffsets)]))
  
  # Primary releases
  
  releases = pdv[, 6:10]
  releases = releases/1000000 #change to MAF
  releases=cbind(rep(NA, nrow(releases)), releases)
  names(releases) = c("Equalization", names(releases[2:ncol(releases)]))
  
  # Create labels
  tier_labs = matrix(NA, nrow = nrow(p.elevation), ncol=6) #hardcoded for now
  tier_labs = as.data.frame(tier_labs)
  
  # Special case for equalization tier
  
  for(k in 1:nrow(tier_labs)){
    
    if(p.elev_delta[k, 1] > 0){
      
      tier_labs[k, 1] = "Equalization"
    }
  }
  
  for(i in 1:nrow(tier_labs)){
    
    for(j in 2:ncol(tier_labs)){
      
      if(is.na(mead_ref[i, j]) == F){
        
        if (minoffsets[i, j] == 0){
          
          tier_labs[i, j] = paste("Release", releases[i, j], "maf;", 
                                  "Balance", releases[i, j] + maxoffsets[i, j], 
                                  "maf if Mead <", mead_ref[i, j], "ft",  sep = " ")
          
        } else{
          
          tier_labs[i, j] = paste("Release", releases[i, j], "maf;", 
                                  "Balance", releases[i, j] + maxoffsets[i, j] - minoffsets[i, j], "-", 
                                  releases[i, j] + maxoffsets[i, j], 
                                  "maf if Mead <", mead_ref[i, j], "ft",  sep = " ")
        }
      }
      
      else if(is.na(releases[i, j])==F) {
        
        tier_labs[i, j] = paste("Release", releases[i,j], "maf;")
        
      }
      
      else{
        
        tier_labs[i, j]=NA
      }
    }
  }
  
  names(tier_labs) = c("Equalization Label", "Tier1 Label", "Tier2 Label", 
                       "Tier3 Label", "Tier4 Label", "Tier5 Label" )
  
  
  # Add dead pool and policy number to every df
  p.elevation$policy = as.character(1:nrow(pdv))
  
  p.elev_delta$dead.pool = 3370
  p.elev_delta$policy = as.character(1:nrow(pdv))
  
  tier_labs$dead.pool = NA
  tier_labs$policy = as.character(1:nrow(pdv))
  
  mead_ref$dead.pool = 3370
  mead_ref$policy = as.character(1:nrow(pdv))
  
  minoffsets$dead.pool = 3370
  minoffsets$policy = as.character(1:nrow(pdv))
  
  maxoffsets$dead.pool = 3370
  maxoffsets$policy = as.character(1:nrow(pdv))
  
  releases$dead.pool = 3370
  releases$policy = as.character(1:nrow(pdv))
  
  # Pivot longer
  
  p.elevation = pivot_longer(p.elevation, cols = 1:ncol(p.elevation) - 1, names_to = 'Tier')
  p.elev_delta = pivot_longer(p.elev_delta, cols = 1:ncol(p.elev_delta) - 1, names_to = 'Tier')
  tier_labs = pivot_longer(tier_labs, cols = 1:ncol(tier_labs) - 1, names_to = 'Tier')
  mead_ref = pivot_longer(mead_ref, cols = 1:ncol(mead_ref) - 1, names_to = 'Tier')
  minoffsets = pivot_longer(minoffsets, cols = 1:ncol(minoffsets) - 1, names_to = 'Tier')
  maxoffsets = pivot_longer(maxoffsets, cols = 1:ncol(maxoffsets) - 1, names_to = 'Tier')
  releases=pivot_longer(releases, cols = 1:ncol(releases) - 1, names_to = 'Tier')
  
  p.df = data.frame(p.elevation, delta = p.elev_delta$value, t_lab = tier_labs$value)
  p.df$Tier = rep(c('a', 'b', 'c', 'd', 'e', 'f', 'g'), length(unique(df$policy)))
  p.df$TierName = rep(c('Equalization', 'Tier1', 'Tier2', 'Tier3', 'Tier4', 'Tier5', 'Dead Pool'), length(unique(p.df$policy)))
  
  # Change to number and add zeros for image ordering in tableau
  p.df$policy = as.numeric(p.df$policy)
  p.df$policy = sprintf("%04d", p.df$policy)
  
  
  ##################### Plotting ##########################
  output_folder <- paste0(output_dir, '//policy_images')
  dir.create(output_folder)
  
  for (i in 1:nrow(condensed_archive)){
    
    policy_id = sprintf("%04d", i)
    
    powell_pol = ggplot(subset(p.df, policy %in% c(policy_id)), 
                        aes(fill = Tier, y = delta, x = policy, label = t_lab)) +  
      
      geom_bar(position = "stack", 
               stat = "identity", 
               color = "black", 
               show.legend = FALSE) +
      
      scale_fill_manual(values=c(
        "#1BBC9B", 
        "#67CC8E",
        "#96ED89",
        "#45BF55",
        "#79BD8F",
        "#289976",
        "#26A69A")) +
      
      geom_text(position = position_stack(vjust = .5), size = 2.5) +
      
      theme_minimal() +
      
      ggtitle("Lake Powell") +
      
      theme(plot.title = element_text(hjust = .5)) +
      
      ylab("PE") +
      
      scale_y_continuous(breaks = seq(3385, 3700, by = 20)) +
      
      coord_cartesian(ylim = c(3385.5, 3700))
    
    
    mead_col = as.character(df$v_col[which(df$policy == policy_id)])
    
    mead_pol = ggplot(data=subset(df, policy %in% c(policy_id)), 
                      aes(fill=Tier, y=delta, x=policy, label=v_lab)) +
      
      geom_bar(position="stack", stat="identity",  color="black",  show.legend = FALSE) +
      
      scale_fill_manual(values = mead_col) +
      
      geom_text(position = position_stack(vjust = .5), size = 2.5) +
      
      theme_minimal() +
      
      ggtitle("Lake Mead") +
      
      theme(plot.title = element_text(hjust = .5)) +
      
      ylab('PE') +
      
      scale_y_continuous(breaks = seq(910, 1220, by = 20)) +
      
      coord_cartesian(ylim=c(910, 1220))
    
    
    powell_pol + mead_pol
    
    both.fig = powell_pol + mead_pol
    
    ggsave(paste0(output_folder, '//', policy_id, ".png"), 
           plot = both.fig, 
           device = "png",
           width = 1650, height = 1420, units = "px",
           dpi = 200)
    
  }
  
  
}



MeadHeatmapMatrix <- function(
    
    condensed_archive, 
    dv_indices,
    max_tiers, 
    max_elev, 
    min_elev, 
    disc_length
    
  ){
  
  # Mead DV indices
  MeadSurplus_idx <- dv_indices$MeadSurplus
  MeadEl_idx <- dv_indices$MeadEl
  MeadV_idx <- dv_indices$MeadV
  
  # For process below to work: if surplus tier does not exist, need to replace 
  # surplus elevation (99999999) with 1200
  condensed_archive[, MeadSurplus_idx][
    condensed_archive[, MeadSurplus_idx] == 99999999
  ] <- 1200
  
  # Initialize heatmap dataframe
  n_policies <- nrow(condensed_archive)
  pool_elevs <- seq(max_elev, min_elev, -1* disc_length)
  pool_elevs <- head(pool_elevs, -1)

  hm_df <- data.frame(matrix(
    nrow = length(pool_elevs),
    ncol = n_policies
  ))
  row.names(hm_df) <- pool_elevs
  
  
  for (i in 1:n_policies){
    
    policy <- condensed_archive[i, ]
    
    # Create elevation and volume vectors to keep track of tiers 
    ## First two entries correspond to surplus and normal tiers
    tier_elevs <- c(max_elev, policy[[MeadSurplus_idx]])
    tier_vols <- c(-1, 0)
    
    ## Fill in remaining tiers
    tier_elevs <- append(tier_elevs, policy[MeadEl_idx])
    tier_vols <- append(tier_vols, policy[MeadV_idx])
    
    tier_elevs <- unname(unlist(tier_elevs))
    tier_vols <- unname(unlist(tier_vols))
    
    # Not all policies have max # of shortage tiers. Check where the last tier
    # occurs.
    if (tail(tier_elevs, n=1) > min_elev){
      idx_last_tier <- length(tier_elevs)
    }else{
      idx_last_tier <- match(min_elev, tier_elevs) - 1
    }
    
    # Delete redundant tiers
    tier_elevs <- c(tier_elevs[1:idx_last_tier], 895)
    tier_vols <- tier_vols[1:idx_last_tier]
    
    # Create vector to store shortage volume at each discretized elevation
    ## -1s for surplus tier and 0s for normal tier.
    vol_at_elev <- c(
      rep(  -1, (tier_elevs[1] - tier_elevs[2]) / disc_length  ),
      rep(   0, (tier_elevs[2] - tier_elevs[3]) / disc_length  )
    )
    
    ## Fill in the associated shortage volumes at each discretized elevation
    ## between bottom of normal tier and deadpool
    for (j in 3:length(tier_vols)){
      vol_at_elev <- append(
        vol_at_elev,
        rep(  tier_vols[j], (tier_elevs[j] - tier_elevs[j + 1]) / disc_length  )
      )
    }
    
    hm_df[, i] <- vol_at_elev
    
  }
  
  return(as.matrix(hm_df))
  
}


SortByStartingShortageElevation <- function(
    
  mead_heatmap_matrix, 
  max_shortage_elev
  
){
  
  # Save original row names
  row_names <- row.names(mead_heatmap_matrix)
  
  # Transpose to use dplyr's arrange function
  hm_df_t <- data.frame(t(mead_heatmap_matrix))
  
  # Sort by each column (elevation) hierarchically
  columns <- names(hm_df_t)
  starting_col <- grep(max_shortage_elev, columns)
  columns_to_sort <- columns[starting_col:length(columns)]
  
  hm_df_t <- hm_df_t %>%
    arrange(across(all_of(columns_to_sort), list(desc)))
  
  # Set original row names, transpose back
  names(hm_df_t) <- row_names
  hm_df <- t(hm_df_t)
  
  return(as.matrix(hm_df))
}


GetMeadColorScheme <- function(){
  shortage_volumes <- c(-1, 0, seq(50, 7000, 10))
  
  reclamation_color_scale <- read.csv("data/vol_gradient.csv")
  color_string <- c("#6d46a5", "#4659a5")
  for (i in 1:nrow(reclamation_color_scale) - 1) {
    color_string <- append(color_string, rep(reclamation_color_scale$color[i], 5))
  }
  color_string <- append(color_string, tail(reclamation_color_scale$color, 1))
  
  # Combine the two lists to create a color scheme/ramp
  colors <- colorRamp2(shortage_volumes, color_string)
  
  return(colors)
  
}



column_mapping <- function(original_col_name){
  
  original_col_name <- modify_string(original_col_name)
  
  col_name = switch(
    original_col_name,
    "Powell_Tier_Elevation_DV.Row.0" = "PT1e",
    "Powell_Tier_Elevation_DV.Row.1" = "PT2e",
    "Powell_Tier_Elevation_DV.Row.2" = "PT3e",
    "Powell_Tier_Elevation_DV.Row.3" = "PT4e",
    "Powell_Tier_Elevation_DV.Row.4" = "PT5e",
    "Powell_Primary_Release_Volume_DV.Row.0" = "PT1Rel",
    "Powell_Primary_Release_Volume_DV.Row.1" = "PT2Rel",
    "Powell_Primary_Release_Volume_DV.Row.2" = "PT3Rel",  
    "Powell_Primary_Release_Volume_DV.Row.3" = "PT4Rel",
    "Powell_Primary_Release_Volume_DV.Row.4" = "PT5Rel",  
    "Powell_Balance_Max_Offset_DV.Row.0" = "MaxOffset1",
    "Powell_Balance_Max_Offset_DV.Row.1" = "MaxOffset2",  
    "Powell_Balance_Max_Offset_DV.Row.2" = "MaxOffset3",
    "Powell_Balance_Max_Offset_DV.Row.3" = "MaxOffset4",  
    "Powell_Balance_Max_Offset_DV.Row.4" = "MaxOffset5",
    "Powell_Balance_Min_Offset_DV.Row.0" = "MinOffset1", 
    "Powell_Balance_Min_Offset_DV.Row.1" = "MinOffset2", 
    "Powell_Balance_Min_Offset_DV.Row.2" = "MinOffset3",  
    "Powell_Balance_Min_Offset_DV.Row.3" = "MinOffset4", 
    "Powell_Balance_Min_Offset_DV.Row.4" = "MinOffset5",  
    "Powell_Mead_Reference_Elevation_DV.Row.0" = "MeadRef1",
    "Powell_Mead_Reference_Elevation_DV.Row.1" = "MeadRef2",            
    "Powell_Mead_Reference_Elevation_DV.Row.2" = "MeadRef3",
    "Powell_Mead_Reference_Elevation_DV.Row.3" = "MeadRef4",
    "Powell_Mead_Reference_Elevation_DV.Row.4" = "MeadRef5",
    "Mead_Surplus_DV" = "surplus_elev",                                   
    "Mead_Shortage_e_DV.Row.0" = "T1e",
    "Mead_Shortage_e_DV.Row.1" = "T2e",                            
    "Mead_Shortage_e_DV.Row.2" = "T3e",
    "Mead_Shortage_e_DV.Row.3" = "T4e",                   
    "Mead_Shortage_e_DV.Row.4" = "T5e",
    "Mead_Shortage_e_DV.Row.5" = "T6e",                   
    "Mead_Shortage_e_DV.Row.6" = "T7e",
    "Mead_Shortage_e_DV.Row.7" = "T8e",                   
    "Mead_Shortage_V_DV.Row.0" = "T1V",
    "Mead_Shortage_V_DV.Row.1" = "T2V",                   
    "Mead_Shortage_V_DV.Row.2" = "T3V",
    "Mead_Shortage_V_DV.Row.3" = "T4V",                   
    "Mead_Shortage_V_DV.Row.4" = "T5V",
    "Mead_Shortage_V_DV.Row.5" = "T6V",                   
    "Mead_Shortage_V_DV.Row.6" = "T7V",
    "Mead_Shortage_V_DV.Row.7" = "T8V",                   
    "Objectives.Avg_Powell_PE" = "Avg.Powell.PE",
    "Objectives.Powell_Release_LTEMP" = "Powell.Release.LTEMP",                   
    "Objectives.Avg_Mead_PE" = "Avg.Mead.PE",
    "Objectives.LB_Shortage_Volume" = "Avg.LB.Shortage",                  
    "Objectives.Avg_Hydrologic_Shortage" = "Avg.Hydrologic.Shortage",
    "Run.Succeeded" = "Run.Succeeded",                               
    "Objectives.LB_Shortage_Frequency" = "LB.Shortage.Frequency",
    "Objectives.Mead_1020" = "Mead.1020",                              
    "Objectives.Powell_3525" = "Powell.3525",
    "Objectives.Powell_3490" = "Powell.3490",                       
    "Objectives.Powell_WY_Release" = "Powell.WY.Release",
    "Objectives.Start_in_EQ" = "Start.In.EQ",                        
    "Objectives.VariableHydrologicShortageIndicatorMetric" = "VariableHydrologicShortageIndicatorMetric",
    "Objectives.Lee_Ferry_Deficit" = "LF.Deficit",                       
    "Objectives.Max_Delta_Annual_Shortage" = "Max.Delta.Shortage",
    "Objectives.Avg_Annual_LB_Policy_Shortage" = "Avg.Policy.Shortage",   
    "Objectives.Mead_1000" = "Mead.1000"                        
  )
  
  if (is.null(col_name)){
    original_col_name
  } else{
    col_name
  }
  
}


modify_string <- function(x) {
  if (grepl("^Objectives.*\\.(1|2)$", x)) {
    return(sub("\\.(1|2)$", "", x))
  } else {
    return(x)
  }
}


get_xml_info <- function(xml_file_path){
  
  xml_file <- read_xml(xml_file_path)
  
  # Get objectives, metrics, and constraints
  objectives <- xml_find_all(xml_file, "//objectiveList/objective/name")
  objective_names <- xml_text(objectives)
  
  constraints <- xml_find_all(xml_file, "//constraintList/constraint/name")
  constraint_names <- xml_text(constraints)
  
  metrics <- xml_find_all(xml_file, "//metricList/metric/name")
  metric_names <- xml_text(metrics)
  
  # Return
  list(objectives = objective_names, constraints = constraint_names, metrics = metric_names)
  
}

get_dv_indices <- function(condensed, col_names, n_powell_tiers, n_mead_tiers){
  
  if (condensed){
    
    PTierEl_first <- which(col_names == "PT1e")
    PRels_first <- which(col_names == "PT1Rel")
    MeadRefs_first <- which(col_names == "MeadRef1")
    BalMaxOffset_first <- which(col_names == "MaxOffset1")
    BalMinOffset_first <- which(col_names == "MinOffset1")
    MeadSurplus <- which(col_names == "surplus_elev")
    MeadEl_first <- which(col_names == "T1e")
    MeadV_first <- which(col_names == "T1V")
    
  } else{
    
    PTierEl_first <- which(col_names == "Powell_Tier_Elevation_DV.Row.0")
    PRels_first <- which(col_names == "Powell_Primary_Release_Volume_DV.Row.0")
    MeadRefs_first <- which(col_names == "Powell_Mead_Reference_Elevation_DV.Row.0")
    BalMaxOffset_first <- which(col_names == "Powell_Balance_Max_Offset_DV.Row.0")
    BalMinOffset_first <- which(col_names == "Powell_Balance_Min_Offset_DV.Row.0")
    MeadSurplus <- which(col_names == "Mead_Surplus_DV")
    MeadEl_first <- which(col_names == "Mead_Shortage_e_DV.Row.0")
    MeadV_first <- which(col_names == "Mead_Shortage_V_DV.Row.0")
    
    
  }
  
  PTierEl <-      c(  PTierEl_first:(PTierEl_first + n_powell_tiers - 1)  )
  PRels <-        c(  PRels_first:(PRels_first + n_powell_tiers - 1)  )
  MeadRefs <-     c(  MeadRefs_first:(MeadRefs_first + n_powell_tiers - 1)  )
  BalMaxOffset <- c(  BalMaxOffset_first:(BalMaxOffset_first + n_powell_tiers - 1)  )
  BalMinOffset <- c(  BalMinOffset_first:(BalMinOffset_first + n_powell_tiers - 1)  )
  MeadEl <-       c(  MeadEl_first:(MeadEl_first + n_mead_tiers - 1)  )
  MeadV <-        c(  MeadV_first:(MeadV_first + n_mead_tiers - 1)  )
  
  # Return
  list(
    PTierEl = PTierEl,
    PRels = PRels,
    MeadRefs = MeadRefs,
    BalMaxOffset = BalMaxOffset,
    BalMinOffset = BalMinOffset,
    MeadSurplus = MeadSurplus,
    MeadEl = MeadEl,
    MeadV = MeadV
  )
  
}





