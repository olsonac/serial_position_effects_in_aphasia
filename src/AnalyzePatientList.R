rm(list=ls())
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

palette_values <- c(
  "#000000","#9b6f02", "#B13401", "#6a0b86","#8862F5", "#0038ae", "#175c12",
  "#339701", "#5dee14", "#DEB604", "#Ee6a14")

shape_values <- c(15,19,17,18,20,3,4,8,10,6)

Sys.setenv(R_CONFIG_ACTIVE = "test")  # uncomment to test
#Sys.setenv(R_CONFIG_ACTIVE = "naming")
#Sys.setenv(R_CONFIG_ACTIVE = "default") 
#Sys.setenv(R_CONFIG_ACTIVE = "no_fractional_values") 

config <- config::get(file="./src/config.yml")
RootDir <- config$root_dir
opts_knit$set(root.dir = config$root_dir)

RandomSamples <- config$random_samples
DoSimulations <- config$do_simulations

#palette_values = as.vector(do.call(config$pals_palette,args=list()))
#palette_values = palette_values[2:length(palette_values)]

# read functions we will use
source(paste0(RootDir,"/src/function_library/sp_functions.R"))

# read list of patients to analyze
ppt_parms <- read.csv(config$patient_param_file)
ppt_length <- nrow(ppt_parms)
RunLabel <- config$run_label
NoFracPres <- config$remove_fractional_values
if("best_model_index_L1" %in% names(ppt_parms)){
  best_model_index_L1 = ppt_parms$best_model_index_L1
}else{
  best_model_index_L1 = rep(config$best_model_default_index_L1,ppt_length)
}

if("best_model_index_L2" %in% names(ppt_parms)){
  best_model_index_L2 = ppt_parms$best_model_index_L2
}else{
  best_model_index_L2 = rep(config$best_model_default_index_L2,ppt_length)
}

if("best_model_index_L3" %in% names(ppt_parms)){
  best_model_index_L3 = ppt_parms$best_model_index_L3
}else{
  best_model_index_L3 = rep(config$best_model_default_index_L3,ppt_length)
}

RunLabel <- config$run_label

for(i in seq(1,nrow(ppt_parms))){
  CurPat <- ppt_parms$patient[i]
  CurTask <- ppt_parms$task[i]
  MinLength <- ppt_parms$min_length[i]
  MaxLength <- ppt_parms$max_length[i]
  BestModelIndexL1 <- best_model_index_L1[i] # in case we need to choose a simpler model close to AIC best
  BestModelIndexL2 <- best_model_index_L2[i]
  BestModelIndexL3 <- best_model_index_L3[i]
  OutputFilename <- paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/reports/",CurPat,"_",CurTask,"_",RunLabel,"_analysis_report.pdf")
  
  RunDir <- paste0(RootDir,"/output/",RunLabel)
  if(!dir.exists(RunDir)){
    dir.create(RunDir,showWarnings = TRUE)
    print(paste0("created run directory: ",RunDir))
  }  
  PatientDir <- paste0(RootDir,"/output/",RunLabel,"/",CurPat)
  if(!dir.exists(PatientDir)){
    dir.create(PatientDir,showWarnings = TRUE)
    print(paste0("created participant directory: ",PatientDir))
  }
  TablesDir <- paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/tables/")
  if(!dir.exists(TablesDir)){
    dir.create(TablesDir,showWarnings = TRUE)
    print(paste0("created tables directory: ",TablesDir))
  }
  FigDir <- paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/fig/")
  if(!dir.exists(FigDir)){
    dir.create(FigDir,showWarnings = TRUE)
    print(paste0("created figures directory: ",FigDir))
  }
  ReportsDir <- paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/reports/")
  if(!dir.exists(ReportsDir)){
    dir.create(ReportsDir,showWarnings = TRUE)
    print(paste0("created reports directory: ",ReportsDir))
  }
  
  # if(!dir.exists(paste0(RootDir,"/output/",RunLabel))){
  #   dir.create(paste0(RootDir,"/output/",RunLabel),showWarnings = TRUE)
  # }  
  # if(!dir.exists(paste0(RootDir,"/output/",RunLabel,"/",CurPat))){
  #   dir.create(paste0(RootDir,"/output/",RunLabel,"/",CurPat),showWarnings = TRUE)
  # }
  # if(!dir.exists(paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/tables/"))){
  #   dir.create(paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/tables/"),showWarnings = TRUE)
  # }
  # if(!dir.exists(paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/fig/"))){
  #   dir.create(paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/fig/"),showWarnings = TRUE)
  # }
  # if(!dir.exists(paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/reports/"))){
  #   dir.create(paste0(RootDir,"/output/",RunLabel,"/",CurPat,"/reports/"),showWarnings = TRUE)
  # }
  
  ReportTitle <- paste(CurPat," - ",CurTask," - ",RunLabel," - Serial position analysis")
  
  # call analysis script
  rmarkdown::render(paste0(RootDir,"/src/AnalyzeOnePatient6_template.Rmd"), 
          params = list(
            CurPat = CurPat,
            CurTask = CurTask,
            MinLength = MinLength,
            MaxLength = MaxLength,
            BestModelIndexL1 = BestModelIndexL1,
            BestModelIndexL2 = BestModelIndexL2,
            BestModelIndexL3 = BestModelIndexL3,
            ReportTitle = ReportTitle,
            RandomSamples = RandomSamples,
            DoSimulations = DoSimulations,
            RemoveFracPres = NoFracPres
  ), output_file = OutputFilename)
}
