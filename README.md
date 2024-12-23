**These data and code are made available under the Open Database License: [http://opendatacommons.org/licenses/odbl/1.0/](https://opendatacommons.org/licenses/odbl/1.0/). Any rights in individual contents of the database are licensed under the Database Contents License: [http://opendatacommons.org/licenses/dbcl/1.0/](https://opendatacommons.org/licenses/dbcl/1.0/)**


###  Quick start: Duplicating the analysis

Further details about the organisation of the repository and running analyses is included below.  

Results contributing to each figure/table can be found in the archive in the locations listed in `[root_dir]/doc/result_location_list.md
`
##### Clone the repository 

- Open a Terminal/Command shell
- Change to the directory where you want the repository located
- Type: `git clone https://github.com/olsonac/serial_position_effects_in_aphasia.git`
- The repository will be in the folder: `serial_position_effects_in_aphasia`
- edit `[root_directory]/src/config.yml` and set the root directory

#####  By position analyses - Analyses of repetition or naming comparing the influence of position, previous errors, previous correct and word length

These are the analyses that form the principle set of results

- start R studio
- Click on `File > Open project` and choose `serial_position_analysis.Rproj`
- Edit the section of `AnalyzePatientList.R` with the `Sys.sentenv()` functions:

```
#Sys.setenv(R_CONFIG_ACTIVE = "test")  # uncomment to test
#Sys.setenv(R_CONFIG_ACTIVE = "naming")
Sys.setenv(R_CONFIG_ACTIVE = "default") 
```

   - Select the line with 'test' for single/small sets of participants (edit `one_patient_parm_file.csv` to add participants and set task to repetition or naming)
   - Select the line 'naming' to analyze naming results
   - Select the line 'default' to analyze repetition results

- Click on `File > Open File` and choose `src/AnalyzePatientList.R`
- Select the entire script by clicking on `Edit` > `Select All`
- Click on `Run` in the bar just above the script

##### By word analyses - Influence of length and frequency

- start R studio
- Click on `File > Open project` and choose `serial_position_analysis.Rproj`
- Click on `File > Open File` and choose `src/ByWord_LenFreqModels.R`
- Select the entire script by clicking on `Edit` > `Select All`
- Click on `Run` in the bar just above the script

##### Plot length x position results for all participants in a grid

- start R studio
- Click on `File > Open project` and choose `serial_position_analysis.Rproj`
- Click on `File > Open File` and choose `src/plot_all_length_pos_in_grid.R`
- Select the entire script by clicking on `Edit` > `Select All`
- Click on `Run` in the bar just above the script

##### Simulate results that depend on position only, previous correct and previous errors.

- start R studio
- Click on `File > Open project` and choose `serial_position_analysis.Rproj`
- Click on `File > Open File` and choose `src/simulation/sim_sp_models.Rmd`
- In the last section, `Test models,` adjust by hand to choose the simulated dataset to examine for interactions (these will occur by chance because simulations always have only one principle data generating factor: word length, previous correct, previous error).  E.g. if examining the simulated data and model fits when the generating function is based on previous correct, choose previous correct for the tables and plots here.  The main (and only) generating factor should be clear, but best models also often include some minor subordinate influences from other factors due to overfitting and this process allows us to see that that is the case.
- In the menu bar just above the script, choose `Run` > `Run All`

# Organization of the repository

There are four directories in the repository:

1) **[root_dir]/data:** Raw data organised in folders for each participant
	- Each participant has one file with position-level data (one row for each position in each stimulus word; used for analysis of position) and one file with word-level data (one row for each word; used for correct/error analyses).  Where there are data from both naming and repetition these are in different files. <br>
	- 
1) **[root_dir]/doc:** Documentation, currently just a copy of this README file <br>

2) **[root_dir]/output:** Results for each participant organized into folders and .csv files for summaries <br>
	- Within each output folder there is:
		- A folder for figures (fig).
		- A folder with the full R output in a .pdf file (reports).
		- A folder with .csv files that have results in tabular format (tables)<br>

1) **[root_dir]/src:** R scripts and .csv files used for control <br>
	- The organization of the analysis scripts is described below.

#### Organziation of the analysis scripts

The core analysis script is ./src/AnalyzeOnePatient6_template.Rmd.  This does the analysis of a single participant.
- The script is called for each patient by a loop in ./src/AnalyzePatientList.R. 

#### Configuring the analysis with `config.yml`

A number of parameters that control the analysis can be set in the file `./src/config.yml`

In the `default:` section:
- `root_dir:`  The top level folder where the repository has been cloned
- `patient_param_file:` The .csv file with the list of participants and analysis parameters for each.
- `best_model_default_index_L1:` The default index for the best model, where models are ordered by AIC value (minimum first).  Normally this is set to 1.  Level 1 is includes just a single explanatory variable.
- `best_model_default_index_L2:` As above.  Level 2 includes models with two explanatory variables.
- `best_model_default_index_L3:` As above,  Level 3 includes models with three explanatory variables
- `random_samples:` The number of permutations to use for permutation testing of cumulative error and cumulative preserved values (currently 100)
- `do_simulations:` TRUE/FALSE value.  If TRUE do permutation testing of cumulative error and cumulative preserved models (simulations take more time, so there are times when this should be FALSE -- e.g. testing something that doesn't impact the simulations section)

In the `test:` section a different set of the parameters listed above can be set.  One common use is to set a different `patient_param_file:` This is typically used, for example, if a single participant needs to be analyzed without doing all of the other participants.  A file that has just a line for the single participant can be added and then the 'test' option can be set in the file `AnalyzePatientList.R` in the following section:

```
#Sys.setenv(R_CONFIG_ACTIVE = "test")  # uncomment to test
#Sys.setenv(R_CONFIG_ACTIVE = "naming")
Sys.setenv(R_CONFIG_ACTIVE = "default") 
```


In the `naming:` section a different set of parameters can be set.  This is typically used to choose a different `patient_parm_file` including just those participants who did the naming task.

In the `patient_parm_file:` (e.g. repetition_patient_list_csv), the following parameters are set:
- `patient` Participant identifier (name of the directory under the data directory that has raw data for this participant, currently a two-letter identifier for each participant)
- `min_length` Sets the minimum word length to analyze.
- `max_length` Sets the maximum word length to analyze.
- `task` (e.g. repetition or naming) - Sets part of the input datafile name (to choose the correct datafile for the task) and part of output specifications (e.g. filenames, the report title)
- [optional] `best_model_index_L2` Sets the index of the model to choose at level 2, where there are two explanatory variables  (1 is the best model. 2 could be a competitive model that is close to 1 and simpler, for example) [This was not needed at level 1 or 3]

#### Testing a limited number of participants instead of the whole set

To test one (or a limited number of) participants, use the 'test' option in AnalyzePatientList.R

i.e. 
```
Sys.setenv(R_CONFIG_ACTIVE = "test")
```

e.g. Currently, wtih this option participant information is taken from the file
```
./src/one_patient_parm_file.csv 
```
(or the file listed in the `patient_parm_file` line of the `test` section of config.yml)

To test a full set of participants, use the 'default' config option:
```
Sys.setenv(R_CONFIG_ACTIVE = "default") 
```
and list all participants in the file 'repetition_patient_list.csv'

#### Functions used by the analysis routines

Functions used by the analysis routines are found in the following file
```
./src/function_library/sp_functions.R
```

This should be loaded by a 'source' command, i.e.:
source(paste0(RootDir,"/src/function_library/sp_functions.R"))

### Running the analysis

To run the analysis, set the config.yml options as needed, adjust which section of the config.yml file is used with `Sys.setenv(R_CONFIG_ACTIVE = "[section name]")`, and run the script AnalyzePatientList.R.  AnalyzePatientList.R will call AnalyzeOnePatient6_template.Rmd once for each participant in the file specified by `patient_parm_file` in config.yml
### Organization of results
output is in `[root_dir]/output`
Each patient has a folder with results. Results are divided into fig/reports/tables
- fig includes plots
- tables includes .csv tables of results
- reports includes the R-markdown report generated by knitting the AnalyzeOnePatient6_template.Rmd for one patient

#### Simulations of serial position effects

The script to make simulation data is in ./src/simulation/sim_sp_models.Rmd
Simulation results are in the same folder (NOT in output)

#### Length X frequency analysis _by word_ (not position)

The script for the length X freq analysis BY WORD is in ./src/ByWord_LenFreqModels.R
This writes model results to each ppt's `[root_dir]/output/tables` folder (with coefficients) and writes a summary file to `[root_dir]/output`

#### Length x position plots for all participants
The script for plotting length x pos results for all participants in a grid is in `[root_dir]/src/plot_all_length_pos_in_grid.R`  

The resulting .tif files are in the main folder:
`[root_dir]/phonological_patients_length_position.tif`
`[root_dir]/mixed_patients_length_position.tif`
`[root_dir]/apraxic_patients_length_position.tif`

#### R Libraries used by `AnalyzePatientList.R`

The following libraries are used by `AnalyzePatientList.R`

```
library(rmarkdown)
library(fmsb) # for TestModels
library(lme4) # for mixed models
library(kableExtra) # for formatting
library(MASS) # for dropterm
library(tidyverse) # for data manipulation and plotting
library(dominanceanalysis) # for dominance analysis
library(cowplot) # for multiple plots in one figure
library(formatters) # to wrap strings for plot titles
library(knitr)
library(pals)
library(ggtext) # for rendering superscripts in plot titles
```


