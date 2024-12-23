library(MASS)
library(tidyverse)
library(lme4)


Sys.setenv(R_CONFIG_ACTIVE = "default") 
config <- config::get(file="./src/config.yml")
RootDir <- config$root_dir


########################################################################
# functions


ModelLenFreqEffects<-function(CurPat){
  WordDatFilename<-paste0("./data/",CurPat,"/",CurPat,"_repetition_worddat.csv")
  cat("reading ",WordDatFilename,"....\n")
  MyData<-read.csv(WordDatFilename)
  ModelSet<-c("correct ~ plen*log_freq","correct ~ plen+log_freq","correct ~ plen","correct ~ log_freq")
  ModelLabels<-c("Interaction","2Main","Len","Freq")
  ModelFamily<-binomial()
  ModelSetResults<-DoGLMModelSet(MyData,ModelSet,ModelLabels,ModelFamily)
  LenImportance<-DoFactorImportance(ModelSetResults,c(1,2,3),"Length")
  FreqImportance<-DoFactorImportance(ModelSetResults,c(1,2,4),"Frequency")
  LFResult<-list(ModelResults=ModelSetResults,LenImportance=LenImportance,FreqImportance=FreqImportance)
  return(LFResult)
}

ExtractCoefficients<-function(ModelResult){
  NumModels <- length(ModelResult)
  coefficient_list <- NULL
  model_formulas <- character(length=length(ModelResult))
  for(Index in seq(1,NumModels)){
    coefficient_list <- c(coefficient_list,names(coefficients(ModelResult[[Index]])))
    model_formulas[Index] <- deparse1(ModelResult[[Index]]$formula)
  }
  coefficient_list <- unique(coefficient_list)
  coefficient_values <- matrix(data=NA,nrow=NumModels,
                               ncol=length(coefficient_list),
                               dimnames=list(model_formulas,coefficient_list))
  for(Index in seq(1,NumModels)){
    model_coefficients <- coefficients(ModelResult[[Index]])
    coefficient_names <- names(model_coefficients)
    rownames(coefficient_values)[Index] <- deparse1(ModelResult[[Index]]$formula)
    for(coefficient_number in seq(1,length(model_coefficients))){
      coefficient_values[deparse1(ModelResult[[Index]]$formula),coefficient_names[coefficient_number]]<-model_coefficients[coefficient_number]
    }
  }
  return(coefficient_values)
}

DoGLMModelSet<-function(MyData,ModelSet,ModelLabels,ModelFamily){
  NumModels<-length(ModelSet)
  ModelResultList<-vector("list",length=NumModels)
  ModelFormulas<-vector("character")
  for(i in seq(1,NumModels)){
    ModelResultList[[i]]<-glm(as.formula(ModelSet[i]),family=ModelFamily,data=MyData)
    ModelFormulas[i]<-deparse1(ModelResultList[[i]]$formula)
  }
  coefficient_values <- ExtractCoefficients(ModelResultList)
  ModelResults<-data.frame(Model=ModelFormulas,
                           AIC=sapply(ModelResultList,AIC),
                           row.names = ModelFormulas)
  MinAIC<-min(ModelResults$AIC)
  ModelResults$DeltaAIC<-ModelResults$AIC-MinAIC
  WeightDenom<-exp(-0.5*ModelResults$DeltaAIC)
  ModelResults$AWeight<-exp(-0.5*ModelResults$DeltaAIC)/sum(WeightDenom)
  MergedModelResults <- merge(ModelResults,coefficient_values, by='row.names',sort=FALSE)
  return(MergedModelResults)
}

DoFactorImportance<-function(ModelResults,ModelIndicies,FactorLabel){
  Importance<-sum(ModelResults$AWeight[ModelIndicies])
  names(Importance)<-FactorLabel
  return(Importance)
}

# End of functions
#############################################################

PhonPat<-c("AC","DS","GM","MC","MP","RM","TC","VS")
MixedPat<-c("AG","CA","MS","PM")
ApraxicPat<-c("AM","AP","AV","DC","DG","EM","GC","MI","OB","PV","SR")

PatList<-c("AC","DS","GM","MC","MP","RM","TC","VS", # phonological
           "AG","CA","MS","PM", # mixed
           "AM","AP","AV","DC","DG","EM","GC","MI","OB","PV","SR") # apraxic
PatClass<-c("phonological",
            "phonological",
            "phonological",
            "phonological",
            "phonological",
            "phonological",
            "phonological",
            "phonological",
            "mixed",
            "mixed",
            "mixed",
            "mixed",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic",
            "apraxic")

############################################################

# calculate length and frequency models using binomial regression on correct/incorrect

LenFreqResults<-lapply(PatList,ModelLenFreqEffects)

for(Index in seq(1,length(PatList))){
  CurPat<-PatList[Index]
  cat("writing coefficients for ppt: ",CurPat,"\n")  
  CurResult <- LenFreqResults[[Index]]$ModelResults # a dataframe with results from all models
  CurResult <- CurResult %>% arrange(AIC) # sort by AIC value
  CurResult <- subset(CurResult, select = -c(Row.names))
  PptFolder <- paste0("./output/",CurPat)
  TablesFolder <- paste0("./output/",CurPat,"/tables")
  CoefficientFilename <- paste0(TablesFolder,"/",CurPat,"_repetition_byword_len_freq_model_results.csv")
  if(!exists(PptFolder)){
    dir.create(PptFolder,showWarnings = TRUE)
  }
  if(!exists(TablesFolder)){
    dir.create(TablesFolder,showWarnings = TRUE)
  }
  write.csv(CurResult,file = CoefficientFilename, row.names = FALSE)
}

InteractionDelta<-sapply(LenFreqResults,with,ModelResults$DeltaAIC[1])
TwoMainDelta<-sapply(LenFreqResults,with,ModelResults$DeltaAIC[2])
LenDelta<-sapply(LenFreqResults,with,ModelResults$DeltaAIC[3])
LogFreqDelta<-sapply(LenFreqResults,with,ModelResults$DeltaAIC[4])
AWInteraction<-sapply(LenFreqResults,with,ModelResults$AWeight[1])
AW2Main<-sapply(LenFreqResults,with,ModelResults$AWeight[2])
AWLen<-sapply(LenFreqResults,with,ModelResults$AWeight[3])
AWLogFreq<-sapply(LenFreqResults,with,ModelResults$AWeight[4])
LenImport<-sapply(LenFreqResults,with,LenImportance)
LogFreqImport<-sapply(LenFreqResults,with,FreqImportance)
PatientLenFreqResults<-data.frame(Patient=as.character(PatList),
                        PType=as.character(PatClass),
                        Interaction=InteractionDelta,
                        LenFreqMainEffects=TwoMainDelta,
                        Len=LenDelta,
                        LogFreq=LogFreqDelta,
                        AWtInteraction=AWInteraction,
                        AWtLenFreqMain=AW2Main,
                        AWtLen=AWLen,
                        AWtLogFreq=AWLogFreq,
                        LenImportance=LenImport,
                        LogFreqImportance=LogFreqImport)
write.csv(PatientLenFreqResults,paste0(RootDir,"/output/all_patients_byword_len_freq_results.csv"),row.names = FALSE)