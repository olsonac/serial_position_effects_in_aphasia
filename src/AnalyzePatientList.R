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

#Sys.setenv(R_CONFIG_ACTIVE = "test")  # uncomment to test
#Sys.setenv(R_CONFIG_ACTIVE = "naming")
Sys.setenv(R_CONFIG_ACTIVE = "default") 

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

for(i in seq(1,nrow(ppt_parms))){
  CurPat <- ppt_parms$patient[i]
  CurTask <- ppt_parms$task[i]
  MinLength <- ppt_parms$min_length[i]
  MaxLength <- ppt_parms$max_length[i]
  BestModelIndexL1 <- best_model_index_L1[i] # in case we need to choose a simpler model close to AIC best
  BestModelIndexL2 <- best_model_index_L2[i]
  BestModelIndexL3 <- best_model_index_L3[i]
  OutputFilename <- paste0(RootDir,"/output/",CurPat,"/reports/",CurPat,"_",CurTask,"_analysis_report.pdf")
  
  if(!dir.exists(paste0("./output/",CurPat))){
    dir.create(paste0("./output/",CurPat),showWarnings = TRUE)
  }
  if(!dir.exists(paste0("./output/",CurPat,"/tables/"))){
    dir.create(paste0("./output/",CurPat,"/tables/"),showWarnings = TRUE)
  }
  if(!dir.exists(paste0("./output/",CurPat,"/fig/"))){
    dir.create(paste0("./output/",CurPat,"/fig/"),showWarnings = TRUE)
  }
  if(!dir.exists(paste0("./output/",CurPat,"/reports/"))){
    dir.create(paste0("./output/",CurPat,"/reports/"),showWarnings = TRUE)
  }
  
  ReportTitle <- paste(CurPat," - ",CurTask," - Serial position analysis (v5)")
  
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
            DoSimulations = DoSimulations
  ), output_file = OutputFilename)
}
