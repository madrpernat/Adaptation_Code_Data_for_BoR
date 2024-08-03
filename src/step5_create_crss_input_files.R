library(readxl)
library(dplyr)
library(magrittr)

# For each SOW in the 500 SOW ensemble, we need input files for CRSS. This
# includes a) flow files for the 29 inflow points. These files were provided by
# Reclamation, b) a UB demand file, and c) initial conditions files. Each SOW
# already has initial elevations for Powell and Mead, and we are going to use
# CRMMS projections to fill in the other initial conditions.


# Load the 500 SOW ensemble information
five_hundred_sow_info = read.csv('output/python_output/500_sow_info.csv')

# Set up paths to read in crmms data
data_dir <- 'data/eocy2026_crmms_data/'
crmms_workbooks <- c(
  paste0(data_dir, 'CrmmsToCrss_Monthly_ESP80.xlsx'), 
  paste0(data_dir, 'CrmmsToCrss_Monthly_ESP90.xlsx'), 
  paste0(data_dir, 'CrmmsToCrss_Monthly_ESP100.xlsx')
)


# Create a dataframe to store EOCY2026 CRMMS run info (i.e., the initial 
# conditions of our post-2026 CRSS runs)

## Extract slot names
slots <- names(data.frame(read_excel(crmms_workbooks[1])))[-1]

## Initialize IC dataframe
ic_df <- data.frame('Slot'=slots)
rownames(ic_df) <- slots


# Initialize dataframe to keep track of the different workbook/worksheet pairs 
# (i.e., the different CRMMS run/trace pairs)
workbook_sheet_index <- data.frame(matrix(nrow=90, ncol = 2))
names(workbook_sheet_index) <- c('Workbook', 'Sheet')

# Populate the IC and workbook/worksheet dataframes
counter <- 1
for (workbook in crmms_workbooks){
  
  sheet_list <- excel_sheets(workbook)
  
  for (sheet in sheet_list){
    
    temp_df <- data.frame(read_excel(workbook, sheet))
    
    temp_df <- temp_df %>% set_rownames(temp_df[,1]) %>% select(-1)
    
    eocy2026 <- t(temp_df["2026-12-01",])
    
    ic_df <- cbind(ic_df, eocy2026)
    
    workbook_sheet_index[counter, ] <- c(workbook, sheet)
    
    counter <- counter + 1
    
  }
  
}

ic_df <- ic_df %>% set_colnames(0:90) %>% select(-1)


# Identify the best matching ICs for each SOW

## Extract Mead and Powell initial pool elevations of 500 SOW ensemble
five_hundred_sow_ic <- five_hundred_sow_info %>% select(c(mead, powell))


## Extract Mead and Powell EOCY2026 pool elevations from CRMMS runs. Change the
## column names to match those of five_hundred_sow_ic
crmms_ic <- data.frame(t(ic_df)) %>% 
  select(c(Mead.Pool.Elevation, Powell.Pool.Elevation)) %>% 
  set_colnames(names(five_hundred_sow_ic))


## For each SOW, identify the CRMMS run whose EOCY2026 Mead and Powell pool
## elevations are closest to those of the SOW (based on Euclidean distance). 
## Save the "best matching" CRMMS run ID for each.
closest_id <- c()
for (i in 1:nrow(five_hundred_sow_ic)){
  
  ## Get Mead and Powell ICs for SOW i
  ic <- five_hundred_sow_ic[i,]
  
  ## Stores the Euclidean distance between the SOW and each CRMMS run
  distances <- apply(crmms_ic, MARGIN=1, function(row) {
    dist(rbind(ic, row))
  })
  
  ## Save CRMMS run index with smallest distance to SOW
  closest_id <- unname(append(closest_id, which.min(distances)))
  
}


# For each SOW, write out initial condition and demand CRSS input files

output_dir <- 'output/r_output/500_sow_crss_input_files/SystemConditionInput'

for (i in 1:nrow(five_hundred_sow_ic)){
  
  folder <- paste0(output_dir, '/trace', i, '/')
  dir.create(folder, recursive = TRUE)
  idx <- closest_id[i]
  workbook_path <- workbook_sheet_index[idx, 1]
  sheet <- workbook_sheet_index[idx, 2]
  
  # Pool elevation slots
  pool_slots <- c(
    "BlueMesa.Pool.Elevation", 
    "Crystal.Pool.Elevation", 
    "Fontenelle.Pool.Elevation",
    "Havasu.Pool.Elevation", 
    "Mead.Pool.Elevation", 
    "Mohave.Pool.Elevation",
    "MorrowPoint.Pool.Elevation", 
    "Navajo.Pool.Elevation", 
    "Powell.Pool.Elevation",
    "TaylorPark.Pool.Elevation"
  )
  for (slot in pool_slots) {
    elevation <- if (slot == "Mead.Pool.Elevation") {
      five_hundred_sow_ic[i, "mead"]
    } else if (slot == "Powell.Pool.Elevation") {
      five_hundred_sow_ic[i, "powell"]
    } else {
      ic_df[slot, idx]
    }
    writeLines(
      c('data_date: 2026-12-31 24:00', 'units: ft', elevation), 
      paste0(folder, sub("\\..*", "", slot), '.poolelevation')
    )
  }
  
  # Bank storage slots
  bank_slots <- c(
    "FlamingGorge.Bank.Storage", 
    "Mead.Bank.Storage", 
    "Powell.Bank.Storage"
  )
  for (slot in bank_slots) {
    writeLines(
      c('data_date: 2026-12-31 24:00', 'units: acre-ft', ic_df[slot, idx]), 
      paste0(folder, sub("\\..*", "", slot), '.bankstorage')
    )
  }
  
  # Additional slots with multiple timesteps
  ## Reopen the excel workbook/sheet associated with the best matching CRMMS run
  temp_df <- data.frame(read_excel(workbook_path, sheet))
  temp_df <- temp_df %>% set_rownames(temp_df[, 1]) %>% select(-1)
  
  ## 3 previous timesteps for Powell.Outflow
  powell_outflows <- temp_df[
    c("2026-10-01", "2026-11-01", "2026-12-01"), 
    "Powell.Outflow"
  ]
  
  ## 6 previous timesteps for GreenRAboveFlamingGorge.Outflow
  green_outflows <- temp_df[
    c("2026-07-01", "2026-08-01", "2026-09-01", "2026-10-01", "2026-11-01", "2026-12-01"), 
    "GreenRAboveFlamingGorge.Outflow"
  ]
  
  ## 6 previous timesteps for FlamingGorge.Pool.Elevation
  flaminggorge_poolelevations <- temp_df[
    c("2026-07-01", "2026-08-01", "2026-09-01", "2026-10-01", "2026-11-01", "2026-12-01"), 
    "FlamingGorge.Pool.Elevation"
  ]
  
  writeLines(
    c('data_date: 2026-10-31 24:00', 'units: acre-ft/month', powell_outflows[1], 
      powell_outflows[2], powell_outflows[3]), 
    paste0(folder,"Powell.outflow")
  )
  
  writeLines(
    c('data_date: 2026-07-31 24:00', 'units: acre-ft/month', green_outflows[1], 
      green_outflows[2], green_outflows[3], green_outflows[4], 
      green_outflows[5], green_outflows[6]), 
    paste0(folder, "GreenRAboveFlamingGorge.outflow")
  )
  
  writeLines(
    c('data_date: 2026-07-31 24:00', 'units: ft', 
      flaminggorge_poolelevations[1], flaminggorge_poolelevations[2], 
      flaminggorge_poolelevations[3], flaminggorge_poolelevations[4], 
      flaminggorge_poolelevations[5], flaminggorge_poolelevations[6]), 
    paste0(folder, "FlamingGorge.poolelevation")
  )
  
  
  # And finally, the UB demand input file
  writeLines(
    c('units: NONE', five_hundred_sow_info[i, "demand"]), 
    paste0(folder, "UB.demand")
  )

}

# For each SOW, locate flow files and copy over

src_dir <- 'data/crss_flow_files'
output_dir <- 'output/r_output/500_sow_crss_input_files/FlowInput'

for (i in 1:nrow(five_hundred_sow_ic)){

  # Identify folder to copy flow files from
  ensemble <- five_hundred_sow_info$Scenario[i]
  trace_number <- five_hundred_sow_info$TraceNumber[i]
  
  folder_name <- switch(
    ensemble,
    "stress_test" = "Stress Test",
    "cmip5" = "CMIP5",
    "npc_adjusted" = "NPC Adjusted",
    "paleo_drought" = "Paleo Drought Resampled",
    "npc_cmip3" = "NPC_MEKO_CMIP3_FULL"
  )
  
  from_folder <- paste0(src_dir, '/', folder_name, '/trace', trace_number, '/')
  to_folder <- paste0(output_dir, '/trace', i, '/')
  dir.create(to_folder, recursive=TRUE)
  
  files <- list.files(from_folder, full.name=TRUE)
  file.copy(from=files, to=to_folder)
  
}
