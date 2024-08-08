# Adaptation_Code_Data_for_BoR
All code and data necessary to replicate the work done in Phases III and IV.

## Cloning the Repository
Only the source code is provided on GitHub due to storage constraints. To access all files and data associated with this project, please follow these steps:

1. **Download the Data:**
   - Visit [this link](https://drive.filen.io/f/41e26773-47f5-4a56-8dbe-a70f542ece2d#az9ZWc4QcZaveQjgtog53xNx0RjgG1o3) to access the data repository.
   - Use the password from the Phase IV report to unlock the folder.
   - Download the folder to your local machine.
   - Extract the downloaded folders.

2. **Prepare the Repository:**
   - Clone the GitHub repository to your local machine:
     ```sh
     git clone https://github.com/madrpernat/Adaptation_Code_Data_for_BoR.git
     ```
   - Place the extracted folders into the cloned repository directory to match the following structure at the root level:
     ```
     Adaptation_Code_Data_for_BoR
     │   Adaptation_Code_Data_for_BoR.R
     │   environment.yaml
     │   README.md
     │
     ├───borg_directories/
     ├───data/
     ├───output/
     ├───riversmart_study__new_experiment_policies_500_sows/
     └───src/
     ```

## Running the Code
The contents of the fully cloned repository are the result of running all scripts (`step01` through `step14`) in the `src` folder, as well as the associated Borg experiments and RiverSMART study (e.g., the SOW ensembles created in `step09` and written to `output/sow_ensembles_for_optimization` were copied over to the corresponding experiment directories in `borg_directories/`, and the CRSS input files created for the 500 SOW ensemble in `step05` were copied to `riversmart_study__new_experiment_policies_500_sows/Model/Inputs` to perform the policy reevaluations).

There are several components of this study that someone may want to change, such as:
   - The number of SOWs sampled in `step04` (i.e., instead of creating a 500 SOW ensemble, one could create a 200 SOW or 1000 SOW ensemble)
   - The linkage method and/or dendrogram cutoff points in `step07`
   - The number of SOWs in each optimization ensemble and/or the number of ensembles created in `step09`
   - The configurations of the various figures created in `step11` and `step14`

### Running the Scripts

1. **Running R Scripts**
   - Assumes the use of RStudio
   - Open RStudio by opening `Adaptation_Code_Data_for_BoR.R`. This will set the working directory to the project's directory.
   - Navigate to the Files tab in the bottom-right pane of the RStudio interface (alongside other tabs such as Plots, Packages, Help) and navigate to the script you'd like to run.
   - You should be able to make your changes and run the script.

2. **Running Python Scripts**
   - Assumes the use of Anaconda for managing environments
   - Use the Anaconda Prompt to create the environment with the project's `environment.yaml` file:
     ```sh
     conda env create -f path_to_project_directory/environment.yaml --name name_the_environment_whatever_you_want
     ```
   - Multiple options for running the scripts:
      1. **Within the Anaconda Prompt**
         - Activate the environment:
           ```sh
           conda activate the_name_that_you_gave_the_environment
           ```
         - Navigate to project directory
           ```sh
           cd path_to_project_directory
           ```
         - Set the `PYTHONPATH` variable temporarily for the current session:
           ```sh
           set PYTHONPATH=.
           ```
         - Run the script you want to run:
           ```sh
           python path_to_script.py
           ```
      2. **Within PyCharm**
         - Open PyCharm and select "Open" from the welcome screen or `File > Open`. Navigate to and select the cloned repository directory (`Adaptation_Code_Data_for_BoR`).
         - Ensure the project is set up to use the newly created environment:
           - Navigate to `File > Settings > Project: Adaptation_Code_Data_for_BoR > Python Interpreter`.
           - Click `Add Interpreter > Add Local Interpreter...`.
           - Within the `Add Python Interpreter` page, click `Conda Environment`.
           - Toggle the `Use existing environment` option, and then use the dropdown to select the environment you created. Click `OK`.
         - You should now be able to run any of the Python scripts by clicking the green "play" button, either at the bottom of the script itself or in the upper right corner of the PyCharm interface.
