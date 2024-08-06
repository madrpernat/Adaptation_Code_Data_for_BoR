source('src/utils/clear.R')
source('src/utils/policy_analysis_library.R')

results_dir <- 'riversmart_study__new_experiment_policies_500_sows/Scenario'
policy_ids_file <- 'riversmart_study__new_experiment_policies_500_sows/policy_ids.csv'
output_file <- 'output/new_experiments_reevaluation_objectives.csv'

riversmart_policy_ids <- read.csv(file = policy_ids_file)

objectives <- c(
  'Avg_Annual_LB_Policy_Shortage',
  'Avg_Mead_PE',
  'Avg_Powell_PE',
  'LB_Shortage_Volume',
  'Lee_Ferry_Deficit',
  'Max_Delta_Annual_Shortage',
  'Mead_1000',
  'Mead_1020',
  'Powell_3490',
  'Powell_3525',
  'Powell_Release_LTEMP',
  'Powell_WY_Release',
  'Start_in_EQ'
)

sows <- 1:500
n_sows <- length(sows)

dfs <- list()

for (i in 1:nrow(riversmart_policy_ids)){
  
  print(i)
  
  experiment <- riversmart_policy_ids$experiment[i]
  policy <- riversmart_policy_ids$policy[i]
  
  # Locate policy's results folder
  policy_id_string <- sprintf('%03d', i)  # Relative to riversmart study
  policy_dir <- (paste0(results_dir, '/policy', policy_id_string))
  
  # Create df to store the policy's values for each objective for each SOW
  policy_df <- data.frame(matrix(
    nrow = n_sows * length(objectives),
    ncol = 5)
  )
  colnames(policy_df) <- c('Experiment', 'Policy', 'SOW', 'Objective', 'Value')
  
  for (objective in objectives){
    
    values <- read_xlsx(
      path = paste0(policy_dir, "/Objectives.", objective, "_Annual.xlsx"), 
      sheet=2
    )
    values <- as.double(values[360,])[-1]  # End of run in row 360
    
    start_idx <- (which(objectives == objective) - 1) * n_sows + 1
    end_idx <- start_idx + n_sows - 1
    
    policy_df[start_idx:end_idx, ] <- data.frame(
      Experiment = rep(experiment, n_sows),
      Policy = rep(policy, n_sows),
      SOW = sows,
      Objective = rep(objective, n_sows),
      Value = values
    )
    
  }
  
  dfs[[i]] <- policy_df
  
}

combined_df <- bind_rows(dfs)

write.csv(
  x = combined_df,
  file = 'output/new_experiments_reevaluation_data.csv',
  row.names = FALSE
)
