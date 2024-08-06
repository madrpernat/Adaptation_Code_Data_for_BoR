source('src/utils/clear.R')
source('src/utils/policy_analysis_library.R')


borg_directories <- c(
  'borg_directories/TieredPowCoordOps_TieredMead_experiment2_seed8_NEW',
  'borg_directories/TieredPowCoordOps_TieredMead_experiment4_seed8_NEW'
)
archive_file_name <- '/borg_run/Archive.txt'
exp_id_file_name <- '/experiment_id.txt'

riversmart_policy_dir <- 
  'riversmart_study__new_experiment_policies_500_sows/Model/Policies/'

# For each policy in the archives, write dv input files for the riversmart study
policy <- c()
experiment <- c()
counter <- 1

for (borg_directory in borg_directories){
  
  archive_df <- read.table(
    file = paste0(borg_directory, archive_file_name),
    header = TRUE
  )
  
  exp_id <- trimws(readLines(
    con = paste0(borg_directory, exp_id_file_name), 
    warn = FALSE
  ))
  
  for (i in 1:nrow(archive_df)){
    
    output_dir <- paste0(riversmart_policy_dir, "policy", counter)
    dir.create(output_dir)
    
    writeLines(c(
      as.character(archive_df$Powell_Tier_Elevation_DV.Row.0[i]),
      as.character(archive_df$Powell_Tier_Elevation_DV.Row.1[i]),
      as.character(archive_df$Powell_Tier_Elevation_DV.Row.2[i]),
      as.character(archive_df$Powell_Tier_Elevation_DV.Row.3[i]),
      as.character(archive_df$Powell_Tier_Elevation_DV.Row.4[i])
      ),
      paste0(output_dir, "/Powell_Tier_Elevation_DV.txt")
    )
    
    writeLines(c(
      as.character(archive_df$Powell_Primary_Release_Volume_DV.Row.0[i]),
      as.character(archive_df$Powell_Primary_Release_Volume_DV.Row.1[i]),
      as.character(archive_df$Powell_Primary_Release_Volume_DV.Row.2[i]),
      as.character(archive_df$Powell_Primary_Release_Volume_DV.Row.3[i]),
      as.character(archive_df$Powell_Primary_Release_Volume_DV.Row.4[i])
      ),
      paste0(output_dir, "/Powell_Primary_Release_Volume_DV.txt")
    )
    
    writeLines(c(
      as.character(archive_df$Powell_Balance_Max_Offset_DV.Row.0[i]),
      as.character(archive_df$Powell_Balance_Max_Offset_DV.Row.1[i]),
      as.character(archive_df$Powell_Balance_Max_Offset_DV.Row.2[i]),
      as.character(archive_df$Powell_Balance_Max_Offset_DV.Row.3[i]),
      as.character(archive_df$Powell_Balance_Max_Offset_DV.Row.4[i])
      ),
      paste0(output_dir, "/Powell_Balance_Max_Offset_DV.txt")
    )
    
    writeLines(c(
      as.character(archive_df$Powell_Balance_Min_Offset_DV.Row.0[i]),
      as.character(archive_df$Powell_Balance_Min_Offset_DV.Row.1[i]),
      as.character(archive_df$Powell_Balance_Min_Offset_DV.Row.2[i]),
      as.character(archive_df$Powell_Balance_Min_Offset_DV.Row.3[i]),
      as.character(archive_df$Powell_Balance_Min_Offset_DV.Row.4[i])
      ),
      paste0(output_dir, "/Powell_Balance_Min_Offset_DV.txt")
    )
    
    writeLines(c(
      as.character(archive_df$Powell_Mead_Reference_Elevation_DV.Row.0[i]),
      as.character(archive_df$Powell_Mead_Reference_Elevation_DV.Row.1[i]),
      as.character(archive_df$Powell_Mead_Reference_Elevation_DV.Row.2[i]),
      as.character(archive_df$Powell_Mead_Reference_Elevation_DV.Row.3[i]),
      as.character(archive_df$Powell_Mead_Reference_Elevation_DV.Row.4[i])
      ),
      paste0(output_dir, "/Powell_Mead_Reference_Elevation_DV.txt")
    )
    
    writeLines(
      as.character(archive_df$Mead_Surplus_DV[i]),
      paste0(output_dir, "/Mead_Surplus_DV.txt")
    )
    
    writeLines(c(
      as.character(archive_df$Mead_Shortage_e_DV.Row.0[i]),
      as.character(archive_df$Mead_Shortage_e_DV.Row.1[i]),
      as.character(archive_df$Mead_Shortage_e_DV.Row.2[i]),
      as.character(archive_df$Mead_Shortage_e_DV.Row.3[i]),
      as.character(archive_df$Mead_Shortage_e_DV.Row.4[i]),
      as.character(archive_df$Mead_Shortage_e_DV.Row.5[i]),
      as.character(archive_df$Mead_Shortage_e_DV.Row.6[i]),
      as.character(archive_df$Mead_Shortage_e_DV.Row.7[i])
      ),
      paste0(output_dir, "/Mead_Shortage_e_DV.txt")
    )
    
    writeLines(c(
      as.character(archive_df$Mead_Shortage_V_DV.Row.0[i]),
      as.character(archive_df$Mead_Shortage_V_DV.Row.1[i]),
      as.character(archive_df$Mead_Shortage_V_DV.Row.2[i]),
      as.character(archive_df$Mead_Shortage_V_DV.Row.3[i]),
      as.character(archive_df$Mead_Shortage_V_DV.Row.4[i]),
      as.character(archive_df$Mead_Shortage_V_DV.Row.5[i]),
      as.character(archive_df$Mead_Shortage_V_DV.Row.6[i]),
      as.character(archive_df$Mead_Shortage_V_DV.Row.7[i])
      ),
      paste0(output_dir, "/Mead_Shortage_V_DV.txt")
    )
    
    counter <- counter + 1
    policy <- append(policy, i)
    experiment <- append(experiment, exp_id)
    
  }
  
  riversmart_policy_ids <- data.frame(
    riversmart_policy_id = 1:length(policy),
    experiment = experiment,
    policy = policy
  )
  
  write.csv(
    x = riversmart_policy_ids,
    file = 'riversmart_study__new_experiment_policies_500_sows/policy_ids.csv',
    row.names = FALSE
  )
  
}
