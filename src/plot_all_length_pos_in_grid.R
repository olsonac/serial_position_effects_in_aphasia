rm(list=ls())
library(tidyverse) # for data manipulation and plotting
library(formatters) # to wrap strings for plot titles
library(ggtext) # for rendering superscripts in plot titles
library(ggpubr) # for ggarrange

palette_values <- c(
  "#000000","#9b6f02", "#B13401", "#6a0b86","#8862F5", "#0038ae", "#175c12",
  "#339701", "#5dee14", "#DEB604", "#Ee6a14")

shape_values <- c(15,19,17,18,20,3,4,8,10,6)

#Sys.setenv(R_CONFIG_ACTIVE = "test")  # uncomment to test
#Sys.setenv(R_CONFIG_ACTIVE = "naming")
Sys.setenv(R_CONFIG_ACTIVE = "default") 

config <- config::get(file="./src/config.yml")
RootDir <- config$root_dir
setwd(RootDir)
# read functions we will use
source(paste0(RootDir,"/src/function_library/sp_functions.R"))

ppt_parms <- read.csv(config$patient_param_file)

num_plots <- nrow(ppt_parms)

plot_list <- list()
plot_list_names <- NULL

for(i in seq(1,num_plots)){
  # set participant and task
  CurPat <- ppt_parms$patient[i]
  CurTask <- ppt_parms$task[i] 
  MinLength <- ppt_parms$min_length[i]
  MaxLength <- ppt_parms$max_length[i]
  
  plot_list_names <- c(plot_list_names,CurPat) # keep a list of ppt names so we
  # can access plots by name (see below where we set names for plot_list)
  
  # load data
  ModelDatFilename<-paste0(RootDir,"/data/",CurPat,"/",CurPat,"_",CurTask,"_posdat.csv")
  if(!file.exists(ModelDatFilename)){
    print("**ERROR** Input data file not found")
    print(sprintf("\tfile: ",ModelDatFilename))
    print(sprintf("\troot dir: ",RootDir))
  }
  PosDat<-read.csv(ModelDatFilename)
  
  PosDat <- PosDat[PosDat$stimlen >= MinLength,] # include only lengths with sufficient N
  PosDat <- PosDat[PosDat$stimlen <= MaxLength,]
  # include only correct and nonword errors (not no response, word or w+nw errors)
  PosDat <- PosDat[(PosDat$err_type == "correct") | (PosDat$err_type == "nonword"), ]
  PosDat <- PosDat[(PosDat$pos) != (PosDat$stimlen), ] 
  # make sure final positions have been removed
  PosDat$stim_number <- NumberStimuli(PosDat)
  
  # do length/position table and to find minimum p(preserved)
  pos_len_summary <- PosDat %>% group_by(stimlen,pos) %>% summarise(preserved = mean(preserved))
  min_preserved_raw <- min(as.numeric(unlist(pos_len_summary$preserved)))
  min_y_value <- 0.6
  if(min_preserved_raw < min_y_value){
    min_y_value <- min_preserved_raw
  }

  obs_linetypes <- c("solid","solid","solid","solid",
                    "solid","solid","solid","solid","solid")
  pred_linetypes <- NULL

  len_pos_plot <- plot_len_pos_obs_predicted(PosDat,
                                             paste0(CurPat), 
                                             NULL, 
                                             c(min_y_value,1.0),
                                             palette_values, 
                                             shape_values,
                                             obs_linetypes, 
                                             pred_linetypes)
  plot_list[[i]]<-len_pos_plot
}
names(plot_list) <- plot_list_names # so we can access plots by ppt name

ncol_plots <- 4
get_nrow_plots <- function(ncol_plots,plot_list_len){
  nrow_plots <- plot_list_len/ncol_plots
  if(nrow_plots > trunc(nrow_plots)){
    nrow_plots <- trunc(nrow_plots) + 1
  }
  return(nrow_plots)
}

# phonological patients
phon_list <- c("AC","DS","GM","MC","MP","RM","TC","VS")
nrow_plots <- get_nrow_plots(ncol_plots,length(phon_list))
phon_plot<-ggarrange(plotlist = plot_list[phon_list],
                     ncol=ncol_plots,nrow=nrow_plots,common.legend = TRUE,legend="right")
ggsave("phonological_patients_length_position.tif",plot=phon_plot,unit="cm",
       width=ncol_plots * (15*.75),height=nrow_plots*(11*.75),compression = "lzw")

# mixed list
mixed_list <- c("AG","CA","MS","PM")
nrow_plots <- get_nrow_plots(ncol_plots,length(mixed_list))
mixed_plot<-ggarrange(plotlist = plot_list[mixed_list],
                     ncol=ncol_plots,nrow=nrow_plots,common.legend = TRUE,legend="right")
ggsave("mixed_patients_length_position.tif",plot=mixed_plot,unit="cm",
       width=ncol_plots * (15*.75),height=nrow_plots*(11*.75),compression = "lzw")

# apraxic patients
apraxic_list <- c("AM","AP","AV","DC","DG","EM","GC","MI","OB","PV","SR")
nrow_plots <- get_nrow_plots(ncol_plots,length(apraxic_list))
apraxic_plot<-ggarrange(plotlist = plot_list[apraxic_list],
                     ncol=ncol_plots,nrow=nrow_plots,common.legend = TRUE,legend="right")
ggsave("apraxic_patients_length_position.tif",plot=apraxic_plot,unit="cm",
       width=ncol_plots * (15*.75),height=nrow_plots*(11*.75),compression = "lzw")

