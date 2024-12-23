CalcCumErrFromPreserved <- function(PosDat){
  AllCumErr <- NULL
  NumWords <- max(PosDat$stim_number)
  for(i in seq(1,NumWords)){
    CurrentPreserved <- PosDat[PosDat$stim_number == i,"preserved"]
    CurrentErrScore <- 1 - CurrentPreserved
    if(length(CurrentErrScore) == 1){
      CurrentErrScore <- 0
    }else{
      CurrentErrScore <- c(0,CurrentErrScore[1:(length(CurrentErrScore)-1)])
    }
    # shift by one because we want PREVIOUS errors, not current errors
    if(length(CurrentPreserved) != length(CurrentErrScore)){
      print(paste0("**ERROR** error on stim number ",i))
    }
    CurrentCumErr <- cumsum(CurrentErrScore)
    AllCumErr <- c(AllCumErr,CurrentCumErr)
  }
  return(AllCumErr)
}

CalcCumPres<-function(PosDat){
  AllCumPres <- NULL
  NumWords <- max(PosDat$stim_number)
  for(i in seq(1,NumWords)){
    CurrentPreserved <- PosDat[PosDat$stim_number == i,"preserved"]
    if(length(CurrentPreserved) == 1){
      CurrentPresScore <- 0
    }else{
      CurrentPresScore <- c(0,CurrentPreserved[1:(length(CurrentPreserved)-1)])
    }
    # shift by one because we want PREVIOUS errors, not current errors
    if(length(CurrentPreserved) != length(CurrentPresScore)){
      print(paste0("**ERROR** error on stim number ",i,
                   " \nCurrentPreserved: ",CurrentPreserved,
                   "\nCurrentPresScore: ",CurrentPresScore))
    }    
    CurrentCumPres <- cumsum(CurrentPresScore)
    AllCumPres <- c(AllCumPres,CurrentCumPres)
  }
  return(AllCumPres)
}

NoLengthOne <- function(PosDat){
  NumWords <- max(PosDat$stim_number)
  GoodItemNumbers <- NULL
  for(i in seq(1,NumWords)){
    CurrentPreserved <- PosDat[PosDat$stim_number == i,"preserved"]
    if(length(CurrentPreserved) > 1){
      GoodItemNumbers <- c(GoodItemNumbers,i)
    }else{
      print(paste0("**WARNING** item number ",i," was length 1"))
    }
  }
  PosDat <- PosDat[PosDat$stim_number %in% GoodItemNumbers,]
  return(PosDat)
}

NumberStimuli <- function(PosDat){
  StimNum <- 0
  NumDataRows<-nrow(PosDat)
  StimNumArray<-numeric(NumDataRows)
  for(i in seq(1,NumDataRows)){
    if(PosDat$pos[i] == 1){
      StimNum <- StimNum + 1
    }
    StimNumArray[i]<-StimNum
  }
  return(StimNumArray)
}

ShuffleWithinWord <- function(PosDat,ColToShuffle){
  AllShuffledValues <- NULL
  NumWords <- max(PosDat$stim_number)
  for(i in seq(1,NumWords)){
    WordDF<-PosDat[PosDat$stim_number == i,]
    ShuffleIndicies <- sample(nrow(WordDF))
    CurShuffledValues <- WordDF[ShuffleIndicies,ColToShuffle]
    AllShuffledValues <- c(AllShuffledValues,CurShuffledValues)
  }
  return(AllShuffledValues)
}

AddReportLine <- function(ReportDF,Description,ParameterValue){
  TempLineDF <- data.frame(Description = Description, ParameterValue = ParameterValue)
  ReportDF <- rbind(ReportDF,TempLineDF)
  return(ReportDF)
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

TestLevel2Models <- function(PosDat, BestFactor, OtherFactor, OtherFactorName=NULL){
  Level2ModelEquations<-c(
    paste0(BestFactor),
    paste0(BestFactor,"+",OtherFactor),
    paste0("preserved ~ ",OtherFactor)
  ) 
  Level2Res<-TestModels(Level2ModelEquations,PosDat)
  BestModelL2<-Level2Res$ModelResult[[ 1 ]]
  BestModelFormulaL2 <- Level2ModelEquations[1]
  
  AICSummary<-data.frame(Model=Level2Res$Model,
                         AIC=Level2Res$AIC,row.names = Level2Res$Model)
  AICSummary$DeltaAIC<-AICSummary$AIC-AICSummary$AIC[1]
  AICSummary$AICexp<-exp(-0.5*AICSummary$DeltaAIC)
  AICSummary$AICwt<-AICSummary$AICexp/sum(AICSummary$AICexp)
  AICSummary$NagR2<-Level2Res$NagR2
  
  AICSummary <- merge(AICSummary,Level2Res$CoefficientValues,
                      by='row.names',sort=FALSE)
  AICSummary <- subset(AICSummary, select = -c(Row.names))
  
  if(is.null(OtherFactorName)){
    OtherFactorName<-OtherFactor
  }
  write.csv(AICSummary,paste0(RootDir,"/output/",CurPat,"/tables/",CurPat,"_",CurTask,
                              "_best_plus_",OtherFactorName,"_models.csv"),row.names = FALSE) 
  return(AICSummary)
}

ConcatAndAddMissingColumns <- function(DF1,DF2){
  Col1 <- names(DF1)
  Col2 <- names(DF2)
  MissingC2 <- setdiff(Col1,Col2)
  MissingC1 <- setdiff(Col2,Col1)
  DF1[MissingC1] <- NA
  DF2[MissingC2] <- NA
  CombDF<-bind_rows(DF1,DF2)
  return(CombDF)
}
  
TestModels<-function(ModelEquations,PosDat,PrintOutput = TRUE, BestLast = FALSE){
  # used by v4
  require(fmsb)
  NumModels<-length(ModelEquations)
  ModelResult<-vector("list",NumModels)
  ModelFormulas<-vector("character")
  AICResult<-vector("numeric")
  ModelLogLikelihood<-vector("numeric")
  RSquare<-vector("numeric")
  NagR2<-vector("numeric")
  CatCor<-vector("numeric")
  PercentCatCor<-vector("numeric")
  DepVar <- trimws(strsplit(ModelEquations[1],"~")[[1]][1])
  NullModel<-glm(as.formula(paste0(DepVar," ~ 1")),family="binomial",data=PosDat)
  NullPredict<-predict(NullModel,type="response")
  NullCatCor<-PredictedCatCorrect(PosDat[DepVar],NullPredict)    
  for(Index in seq(1,NumModels)){
    ModelResult[[Index]]<-glm(as.formula(ModelEquations[Index]),
                              family="binomial",data=PosDat)
    ModelFormulas[Index]<-deparse1(ModelResult[[Index]]$formula)
    # CoefficientValues[Index,]<-ExtractCoefficients(ModelResult[[Index]],CoefficientList)
    
    # AICResult[Index]<-ModelResult[[Index]]$aic
    AICResult[Index]<-AIC(ModelResult[[Index]])
    ModelLogLikelihood[Index]<-ModelResult[[Index]]$deviance/-2
    PredictedData<-predict(ModelResult[[Index]],type="response")
    NagR2Result<-NagelkerkeR2(ModelResult[[Index]])
    NagR2[Index]<-unlist(NagR2Result$R2)
    CatCor[Index]<-PredictedCatCorrect(PosDat[DepVar],PredictedData)
    PercentCatCor[Index]<-(CatCor[Index]-NullCatCor)/(1-NullCatCor)
    # range that is predictable is 1-Null since Null sets lowest possible prediction
    # % calculated here is the improvement of a model over null as a portion
    # of 1-Null.  This is the portion of the possible predictable range that
    # is actually covered by the model.
  }
  coefficient_values <- ExtractCoefficients(ModelResult)
  AICOrder<-order(AICResult,decreasing = BestLast)
  if(PrintOutput){
    for(Index in seq(1,NumModels)){
      cat("***************************\n")
      cat("model index: ",AICOrder[Index],"\n")
      print(ModelResult[[ AICOrder[Index] ]])
      cat("log likelihood: ",ModelLogLikelihood[ AICOrder[Index] ],"\n")
      cat("Nagelkerke R2: ",NagR2[ AICOrder[Index] ],"\n")
      cat("% pres/err predicted correctly: ",CatCor[ AICOrder[Index] ],"\n")
      cat("% of predictable range [ (model-null)/(1-null) ]: ", PercentCatCor[ AICOrder[Index] ],"\n")
    }
    cat("***************************\n")
  }
  ReturnValue<-list(ModelResult=ModelResult[ AICOrder ],
                    Model=ModelFormulas[ AICOrder ],
                    AIC=AICResult[ AICOrder ],
                    AICOrder=AICOrder,
                    LogLikelihood=ModelLogLikelihood[ AICOrder ],
                    CoefficientValues=coefficient_values[ AICOrder, ],
                    NagR2=NagR2[ AICOrder ],
                    PredCor=CatCor[ AICOrder ],
                    PercentPredCor=PercentCatCor[ AICOrder ])
  return(ReturnValue)
}

EvaluateSubsetData <- function(PosDatSubset,ModelEquations){
  SubsetModelRes<-TestModels(ModelEquations,PosDatSubset)
  BestMEModel<-SubsetModelRes$ModelResult[[ 1 ]]
  BestMEModelFormula<-SubsetModelRes$Model[[1]]
  
  SubsetModelAICSummary<-data.frame(Model=SubsetModelRes$Model,
                                   AIC=SubsetModelRes$AIC,
                                   row.names=SubsetModelRes$Model)
  SubsetModelAICSummary$DeltaAIC <- SubsetModelAICSummary$AIC -
    SubsetModelAICSummary$AIC[1]
  SubsetModelAICSummary$AICexp<-exp(-0.5*SubsetModelAICSummary$DeltaAIC)
  SubsetModelAICSummary$AICwt<-SubsetModelAICSummary$AICexp/
    sum(SubsetModelAICSummary$AICexp)
  SubsetModelAICSummary$NagR2<-SubsetModelRes$NagR2
  
  SubsetModelAICSummary <- merge(SubsetModelAICSummary,
                                SubsetModelRes$CoefficientValues,
                                by='row.names',sort=FALSE)
  SubsetModelAICSummary <- subset(SubsetModelAICSummary, select = -c(Row.names))
  return(SubsetModelAICSummary)
}

PredictedCatCorrect<-function(Obs,Pred){
  # using abs diff
  RetValue<-1-(sum(abs(Obs-Pred))/length(Obs))
  return(RetValue)
}

plot_len_pos_obs_predicted <- function(MyPosDat,PlotTitle,
                                       FittedVariableName=NULL,
                                       min_max=NULL, palette_values = NULL, 
                                       shape_values = NULL,
                                       obs_linetypes=NULL,pred_linetypes=NULL){
  # len/pos table
  if(!is.null(FittedVariableName)){
    fitted_pos_len_summary <- MyPosDat %>% group_by(stimlen,pos) %>% summarise(fitted = mean(!!as.symbol(FittedVariableName)))
    fitted_pos_len_summary$stimlen<-factor(fitted_pos_len_summary$stimlen)
    fitted_pos_len_summary$pos<-as.numeric(as.character((fitted_pos_len_summary$pos)))
  }
  pos_len_summary <- MyPosDat %>% group_by(stimlen,pos) %>% summarise(preserved = mean(preserved))
  #  fitted_pos_len_table <- fitted_pos_len_summary %>% pivot_wider(names_from = pos, values_from = fitted)
  pos_len_summary_N <- MyPosDat %>% group_by(stimlen,pos) %>% summarise(N = n())
  pos_len_table_N <- pos_len_summary_N %>% pivot_wider(names_from = pos, values_from = N)
  pos_len_summary$stimlen <- as.factor(pos_len_summary$stimlen)
  pos_len_summary$pos <- as.numeric(as.character(pos_len_summary$pos))
  
  if(is.null(obs_linetypes)){
    obs_linetypes = "solid"
  }
  if(is.null(pred_linetypes)){
    pred_linetypes = "dashed"
  }
  
  fitted_len_pos_plot <- ggplot(pos_len_summary,
                                aes(x=pos,y=preserved,
                                    group=stimlen,
                                    shape=stimlen,
                                    color=stimlen,
                                    linetype=stimlen)) + 
    geom_point() + geom_line(size=1,alpha=0.8)
  if(!is.null(FittedVariableName)){
    fitted_len_pos_plot <- fitted_len_pos_plot +
    geom_line(data=fitted_pos_len_summary,
              aes(x=pos,y=fitted,
                  group=stimlen,
                  color=stimlen),linetype = pred_linetypes,size=1,alpha=0.5) +
    geom_point(data=fitted_pos_len_summary,
               aes(x=pos,y=fitted,
                   group=stimlen,
                   color=stimlen,
                   shape=stimlen))
  }
  
  fitted_len_pos_plot <- fitted_len_pos_plot + ggtitle(PlotTitle) 
  fitted_len_pos_plot <- fitted_len_pos_plot + 
    scale_x_continuous(name="Word position",breaks=c(1,2,3,4,5,6,7,8,9,10,11)) +
    scale_y_continuous(name = "p(preserved)")
    theme(legend.position="none")
  fitted_len_pos_plot <- fitted_len_pos_plot + 
    scale_shape_manual(name="Word length",values = shape_values) +
    scale_linetype_manual(name="Word length",values=obs_linetypes)
  if(!is.null(palette_values)){
      fitted_len_pos_plot <- fitted_len_pos_plot + 
        scale_color_manual(name="Word length",values = palette_values)
  }
  if(!is.null(min_max)){
    fitted_len_pos_plot <- fitted_len_pos_plot + ylim(min_max)
  }
  return(fitted_len_pos_plot)
}

PlotObsPreviousCorrect<-function(MyData,palette_values,shape_values){
  MyPlot<-PlotPreviousCorrect(MyData,"preserved",palette_values,shape_values)
  return(MyPlot)
}

PlotPreviousCorrect<-function(MyData,ColName,palette_values = NULL,shape_values = NULL){
  LabelSize<-16
  MyData$CumPres<-round(MyData$CumPres)  
  PreCorTable<-MyData %>% group_by(CumPres,pos) %>% summarise(Pres=mean(get(ColName)),N=length(get(ColName)))
  PreCorTable$CumPres<-as.factor(PreCorTable$CumPres)
  PreCorTable<-PreCorTable[PreCorTable$N>10,]
  
  CorPlot<-ggplot(PreCorTable,aes(x=pos,y=Pres,group=CumPres,
                                  shape=CumPres,
                                  color=CumPres))+
    geom_point(size=2)+geom_line(size=1,linetype="solid")
  CorPlot <- FormatPCorPErrPlots(CorPlot,"Number\npreviously\ncorrect", 
                                 LabelSize, palette_values, shape_values)
  return(CorPlot)
}

FormatPCorPErrPlots <- function(CurrentPlot,LegendTitle,
                                LabelSize,
                                palette_values=NULL,
                                shape_values=NULL){
  # CurrentPlot<-CurrentPlot + theme_bw()
  # CurrentPlot<-CurrentPlot + scale_linetype_manual(name=LegendTitle, 
  #                                          values=c("solid","dashed","dotdash","dotted","longdash",
  #                                                   "twodash","solid","dashed","dotdash","dotted","longdash"))
  if(!is.null(shape_values)){
    CurrentPlot<-CurrentPlot + scale_shape_manual(name=LegendTitle,values = shape_values)
  }else{
    CurrentPlot<-CurrentPlot + scale_shape_manual(name=LegendTitle,values = c(15,19,17,18,20,3,4,8,10,6))
  }
  if(!is.null(palette_values)){
    CurrentPlot <- CurrentPlot + scale_color_manual(name=LegendTitle, values = palette_values)
  }
  CurrentPlot<-CurrentPlot + scale_x_continuous(name="Word position",breaks=c(1,2,3,4,5,6,7,8,9,10,11))
  CurrentPlot<-CurrentPlot + scale_y_continuous(name="p(preserved)")
  # change axis label and legend label size
  CurrentPlot<-CurrentPlot + theme(axis.title.x = element_text(size=LabelSize))
  CurrentPlot<-CurrentPlot + theme(axis.title.y = element_text(size=LabelSize))
  # CurrentPlot<-CurrentPlot + theme(legend.title = element_text(size=LabelSize))
  CurrentPlot<-CurrentPlot + theme(legend.text = element_text(size=LabelSize-2))
  CurrentPlot<-CurrentPlot + theme(legend.title = element_text(size=LabelSize-2))
  # change axis tick font size
  CurrentPlot<-CurrentPlot + theme(axis.text.x = element_text(size=LabelSize))
  CurrentPlot<-CurrentPlot + theme(axis.text.y = element_text(size=LabelSize))

  return(CurrentPlot)
}

PlotObsPreviousError<-function(MyData,palette_values,shape_values){
  MyPlot<-PlotPreviousError(MyData,"preserved",palette_values,shape_values)
  return(MyPlot)
}

PlotPreviousError<-function(MyData,ColName,palette_values=NULL,shape_values=NULL){
  LabelSize<-16
  MyData$CumErr<-round(MyData$CumErr)  
  PreErrTable<-MyData %>% group_by(CumErr,pos) %>% summarise(Pres=mean(get(ColName)),N=length(get(ColName))) 
  PreErrTable<-PreErrTable[PreErrTable$N>10,]
  PreErrTable$CumErr<-as.factor(PreErrTable$CumErr)
  
  ErrPlot<-ggplot(PreErrTable,aes(x=pos,y=Pres,group=CumErr,
                                  shape=CumErr,
                                  color=CumErr))+
    geom_point(size=2)+geom_line(size=1,linetype="solid")
  ErrPlot <- FormatPCorPErrPlots(ErrPlot,"Number\nof previous\nerrors",
                                 LabelSize,palette_values,shape_values)
  return(ErrPlot)
}

PlotObsPredPreviousCorrect<-function(MyData,Obs,Predicted,palette_values,shape_values){
  MyPlot<-PlotTwoVarPreviousCorrect(MyData,Obs,Predicted,palette_values,palette_values,shape_values)
  return(MyPlot)
}

PlotTwoVarPreviousCorrect<-function(MyData,ColName1,ColName2,
                                    color1=palette_values,
                                    color2=palette_values,
                                    shape_values){
  LabelSize<-16
  MyData$CumPres<-round(MyData$CumPres)  
  PreCorTableV1<-MyData %>% group_by(CumPres,pos) %>% summarise(Pres=mean(get(ColName1)),N=length(get(ColName1)))
  PreCorTableV1$CumPres<-as.factor(PreCorTableV1$CumPres)
  PreCorTableV1<-PreCorTableV1[PreCorTableV1$N>10,]
  
  CorPlot<-ggplot(PreCorTableV1,aes(x=pos,y=Pres,
                                    group=CumPres,
                                    shape=CumPres,
                                    color=CumPres))+
    geom_point(size=3)+
    geom_line(size=1,linetype="solid")

  PreCorTableV2<-MyData %>% group_by(CumPres,pos) %>% 
    summarise(Pres=mean(get(ColName2)),N=length(get(ColName2)))
  PreCorTableV2$CumPres<-as.factor(PreCorTableV2$CumPres)
  PreCorTableV2<-PreCorTableV2[PreCorTableV2$N>10,]
  
  CorPlot<-CorPlot + geom_point(data=PreCorTableV2,
                                aes(x=pos,y=Pres,
                                    group=CumPres,
                                    shape=CumPres,
                                    color=CumPres),
                                size=2) +
    geom_line(data=PreCorTableV2,aes(x=pos,y=Pres,
                                     group=CumPres,
                                     color=CumPres),
              size=1,linetype="dashed",alpha=0.4) +
    geom_point(data=PreCorTableV2,aes(x=pos,y=Pres,
                                      group=CumPres,
                                      color=CumPres),size=2)
  CorPlot <- FormatPCorPErrPlots(CorPlot,"Number\npreviously\ncorrect",
                                 LabelSize,palette_values,shape_values)  
  return(CorPlot)
}


PlotObsPredPreviousError<-function(MyData,Obs,Predicted,palette_values,shape_values){
  MyPlot<-PlotTwoVarPreviousError(MyData,Obs,Predicted,palette_values,palette_values,shape_values)
  return(MyPlot)
}

PlotTwoVarPreviousError<-function(MyData,ColName1,ColName2,
                                  color1=palette_values,color2=palette_values,
                                  shape_values){
  LabelSize<-16
  MyData$CumErr<-round(MyData$CumErr)  
  PreErrTableV1<-MyData %>% group_by(CumErr,pos) %>% summarise(Pres=mean(get(ColName1)),N=length(get(ColName1))) 
  PreErrTableV1<-PreErrTableV1[PreErrTableV1$N>10,]
  PreErrTableV1$CumErr<-as.factor(PreErrTableV1$CumErr)
  
  ErrPlot<-ggplot(PreErrTableV1,aes(x=pos,y=Pres,
                                    group=CumErr,
                                    shape=CumErr,
                                    color=CumErr))+
    geom_point(size=2)+geom_line(size=1,linetype="solid")

  PreErrTableV2<-MyData %>% group_by(CumErr,pos) %>% summarise(Pres=mean(get(ColName2)),N=length(get(ColName2))) 
  PreErrTableV2<-PreErrTableV2[PreErrTableV2$N>10,]
  PreErrTableV2$CumErr<-as.factor(PreErrTableV2$CumErr)

  ErrPlot<-ErrPlot + geom_point(data=PreErrTableV2,aes(x=pos,y=Pres,
                                                       group=CumErr,
                                                       shape=CumErr,
                                                       color=CumErr),
                                size=2)+
    geom_line(data=PreErrTableV2,aes(x=pos,y=Pres,
                                     group=CumErr,
                                     color=CumErr),
              size=1,linetype="dashed",alpha=0.4)  
  
  ErrPlot <- FormatPCorPErrPlots(ErrPlot,"Number\nof previous\nerrors",
                                 LabelSize,palette_values,shape_values)  
  return(ErrPlot)
}

ModelSetFromDeletions<-function(BestModelFormula,BestModelDeletions,NumModels){
  ModelPrefix <- paste0(strsplit(BestModelFormula,"~")[[1]][1],"~ ")
  BestModelDeletions <- BestModelDeletions %>% arrange(desc(AIC))
  BestModelDeletions <- BestModelDeletions[row.names(BestModelDeletions) != "<none>",]
  # desc = descending -- 
  # model with biggest AIC difference when term is deleted is where we start
  # this is the factor that makes the biggest difference when it is eliminated
  ModelSet<-character(NumModels) # save space for NumModels formulas
  CurrentModel<-ModelPrefix
  for(i in seq(1,NumModels)){ 
    ModelSet[i] <- paste0(CurrentModel, row.names(BestModelDeletions)[i])
    CurrentModel <- paste0(ModelSet[i],"+")
  }
  return(ModelSet)
}

PlotModelSetWithFreq<-function(PosDat,ModelEquations,
                               palette_values=NULL,
                               ModelTitles=NULL,PptID=NULL){
  # ModelEquations should be in order of larger to smaller influences
  # based on AIC for single deletions from best model
  NumModels<-length(ModelEquations)
  NumPlotRows <- NumModels + 1 # plot first line with observed data only
  for(i in seq(1,length(ModelTitles))){
    # translate model variable names from code into readable items
    ModelTitles[i] <- gsub("preserved ~ ","",ModelTitles[i])
    ModelTitles[i] <- gsub("I\\(pos\\^2\\)","position2",ModelTitles[i])
    ModelTitles[i] <- gsub("CumPres","cumulative correct",ModelTitles[i])
    ModelTitles[i] <- gsub("CumErr","cumulative error",ModelTitles[i])
    ModelTitles[i] <- gsub("log_freq","log(frequency)",ModelTitles[i])
    ModelTitles[i] <- gsub("stimlen","word length",ModelTitles[i])
    ModelTitles[i] <- gsub("\\+"," + ",ModelTitles[i])
    ModelTitles[i] <- wrap_string(ModelTitles[i],25,collapse="\n")
    ModelTitles[i] <- gsub("position2","position<sup>2</sup>",ModelTitles[i])
  }
  ModelTitles <- c("Observed data only",ModelTitles)
  PosDatNumCol<-ncol(PosDat)
  ModelResults<-TestModels(ModelEquations,PosDat,BestLast=TRUE)
  PlotList<-list()
  NumTitlePlots<-4 # this is based on the 
  # number of plot columns (len, cumcor, cumerr, freq)
  # and is the number of plots with titles at the top of the column
  RelHeights<-c(0.25,rep(1,NumPlotRows))
  
  # save title spaces -- these are in the first line of plots
  # one at a time (rather than loop) because we may need to customize
  PlotList[[1]]<-ggdraw()
  PlotList[[2]]<-ggdraw()
  PlotList[[3]]<-ggdraw()
  PlotList[[4]]<-ggdraw()
  
  for(i in seq(1,NumModels+1)){ # plus one because the first line is obs data only
    if(i == 1){
      PlotPredicted = FALSE # first line is observed only 
    }else{
      PlotPredicted = TRUE
    }
    if(i == 2){
      PlotLegend<-FALSE # change this if you want a legend
    }else{
      PlotLegend<-FALSE
    }
    if(i > 1){ # first i is observed only
      PosDat[,PosDatNumCol+i-1]<-fitted(ModelResults$ModelResult[[i-1]])
    }
    PlotNumber<-NumTitlePlots+((i-1)*4)
    # skip first line of plots which is titles only (see above)
    CurrentModelTitle<-ModelTitles[i]
    PlotList[[PlotNumber+1]]<-PlotPos(PosDat,PosDatNumCol+i-1,palette_values,
                                      PlotPredicted,PlotLegend,
                                      PlotTitle=CurrentModelTitle)
    #    if(i == NumModels){
    #      PlotList[[PlotNumber+1]]<-add_sub(PlotList[[PlotNumber+1]],"Length")
    #    }
    PlotList[[PlotNumber+2]]<-PlotCumCor(PosDat,PosDatNumCol+i-1,palette_values,
                                         PlotPredicted,PlotLegend)
    PlotList[[PlotNumber+3]]<-PlotCumErr(PosDat,PosDatNumCol+i-1,palette_values,
                                         PlotPredicted,PlotLegend)
    PlotList[[PlotNumber+4]]<-PlotFreq(PosDat,PosDatNumCol+i-1,palette_values,
                                       PlotPredicted,PlotLegend)
    # add here for more columns and change NumTitlePlots
  }
  OAPlot<-plot_grid(plotlist=PlotList,align="h",labels=NULL,ncol=NumTitlePlots,
                    rel_heights=RelHeights)
  OAPlot<-OAPlot+draw_label("Length",x=0.15,y=.96,size=10)
  OAPlot<-OAPlot+draw_label("Previous correct",x=0.41,y=0.96,size=10)
  OAPlot<-OAPlot+draw_label("Previous error",x=0.66,y=0.96,size=10)
  OAPlot<-OAPlot+draw_label("Frequency",x=0.91,y=0.96,size=10)
  OAPlot<-OAPlot+draw_label(PptID,x=0.48,y=.99,size=10)
  return(OAPlot)
}

PlotModelSetWithFreq<-function(PosDat,ModelEquations,
                               palette_values=NULL,
                               ModelTitles=NULL,PptID=NULL){
  # ModelEquations should be in order of larger to smaller influences
  # based on AIC for single deletions from best model
  NumModels<-length(ModelEquations)
  NumPlotRows <- NumModels + 1 # plot first line with observed data only
  for(i in seq(1,length(ModelTitles))){
    # translate model variable names from code into readable items
    ModelTitles[i] <- gsub("preserved ~ ","",ModelTitles[i])
    ModelTitles[i] <- gsub("I\\(pos\\^2\\)","position2",ModelTitles[i])
    ModelTitles[i] <- gsub("CumPres","cumulative correct",ModelTitles[i])
    # changed to position<sup>2</sup> below after line is split (otherwise <sup> confuses splitting)
    ModelTitles[i] <- gsub("CumErr","cumulative error",ModelTitles[i])
    ModelTitles[i] <- gsub("log_freq","log(frequency)",ModelTitles[i])
    ModelTitles[i] <- gsub("stimlen","word length",ModelTitles[i])
    ModelTitles[i] <- gsub("\\+"," + ",ModelTitles[i])
    ModelTitles[i] <- wrap_string(ModelTitles[i],50,collapse="\n")
    ModelTitles[i] <- gsub("position2","position<sup>2</sup>",ModelTitles[i])
    ModelTitles[i] <- gsub("\n","<br>",ModelTitles[i])
  }
  ModelTitles <- c("Observed data only",ModelTitles)
  PosDatNumCol<-ncol(PosDat)
  ModelResults<-TestModels(ModelEquations,PosDat,BestLast=TRUE)
  PlotList<-list()
  NumTitlePlots<-4 # this is based on the 
  # number of plot columns (len, cumcor, cumerr, freq)
  # and is the number of plots with titles at the top of the column
  RelHeights<-c(0.25,rep(1,NumPlotRows))
  
  # save title spaces -- these are in the first line of plots
  # one at a time (rather than loop) because we may need to customize
  PlotList[[1]]<-ggdraw()
  PlotList[[2]]<-ggdraw()
  PlotList[[3]]<-ggdraw()
  PlotList[[4]]<-ggdraw()
  
  for(i in seq(1,NumModels+1)){ # plus one because the first line is obs data only
    if(i == 1){
      PlotPredicted = FALSE # first line is observed only 
    }else{
      PlotPredicted = TRUE
    }
    if(i == 2){
      PlotLegend<-FALSE # change this if you want a legend
    }else{
      PlotLegend<-FALSE
    }
    if(i > 1){ # first i is observed only
      PosDat[,PosDatNumCol+i-1]<-fitted(ModelResults$ModelResult[[i-1]])
    }
    PlotNumber<-NumTitlePlots+((i-1)*4)
    # skip first line of plots which is titles only (see above)

    PlotList[[PlotNumber+1]]<-PlotPos(PosDat,PosDatNumCol+i-1,palette_values,
                                      PlotPredicted,PlotLegend,
                                      PlotTitle=ModelTitles[i])
    #    if(i == NumModels){
    #      PlotList[[PlotNumber+1]]<-add_sub(PlotList[[PlotNumber+1]],"Length")
    #    }
    PlotList[[PlotNumber+2]]<-PlotCumCor(PosDat,PosDatNumCol+i-1,palette_values,
                                         PlotPredicted,PlotLegend)
    PlotList[[PlotNumber+3]]<-PlotCumErr(PosDat,PosDatNumCol+i-1,palette_values,
                                         PlotPredicted,PlotLegend)
    PlotList[[PlotNumber+4]]<-PlotFreq(PosDat,PosDatNumCol+i-1,palette_values,
                                       PlotPredicted,PlotLegend)
    # add here for more columns and change NumTitlePlots
  }
  OAPlot<-plot_grid(plotlist=PlotList,align="h",labels=NULL,ncol=NumTitlePlots,
                    rel_heights=RelHeights)
  OAPlot<-OAPlot+draw_label("Length",x=0.15,y=.96,size=10)
  OAPlot<-OAPlot+draw_label("Previous correct",x=0.41,y=0.96,size=10)
  OAPlot<-OAPlot+draw_label("Previous error",x=0.66,y=0.96,size=10)
  OAPlot<-OAPlot+draw_label("Frequency",x=0.91,y=0.96,size=10)
  OAPlot<-OAPlot+draw_label(PptID,x=0.48,y=.99,size=10)
  return(OAPlot)
}

PlotPos<-function(PosDat,PredCol,palette_values,
                  PlotPredicted=TRUE,PlotLegend=TRUE,PlotTitle=NULL){
  # observed data
  ObsTablePos<-PosDat %>% group_by(stimlen,pos) %>% summarise(Pres=mean(preserved),N=length(preserved))
  ObsTablePos$stimlen<-as.factor(ObsTablePos$stimlen)
  ObsTablePos<-ObsTablePos[ObsTablePos$N>10,]
  
  if(PlotPredicted){
    PredTablePos<-PosDat %>% group_by(stimlen,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(mean))
    PredTableN<-PosDat %>% group_by(stimlen,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(length))
    names(PredTablePos)[length(PredTablePos[1,])]<-"PredPres"
    names(PredTableN)[length(PredTableN[1,])]<-"N"
    PredTablePos$stimlen<-as.factor(PredTablePos$stimlen)
    PredTablePos<-PredTablePos[PredTableN$N>10,]
  }
  
  # plot observed
  PosPlot<-ggplot(ObsTablePos,aes(x=pos,y=Pres,
                                  group=stimlen,
                                  shape=stimlen,
                                  color=stimlen))+
    geom_point(size=2)+geom_line(size=1,linetype="solid")
  
  # plot predicted
  if(PlotPredicted){
    PosPlot<-PosPlot+geom_point(data=PredTablePos,
                                aes(x=pos,y=PredPres,
                                    group=stimlen,
                                    shape=stimlen,
                                    color=stimlen),
                                size=2)+
      geom_line(data=PredTablePos,
                aes(x=pos,y=PredPres,
                    group=stimlen,
                    color=stimlen),linetype="dashed",size=1.1)
  }
  
  PosPlot<-PosPlot + scale_shape_discrete(name="Word length")
  # PosPlot<-PosPlot + scale_linetype_discrete(name="Word length")
  PosPlot<-PosPlot + scale_x_continuous(name="Word position",
                                        breaks=c(1,2,3,4,5,6,7,8,9,10))  +
    theme(axis.title.x = element_text(size=8))
  PosPlot<-PosPlot + scale_y_continuous(name="p(preserved)") +
    theme(axis.title.y = element_text(size=8))
  PosPlot <- PosPlot + scale_color_manual(values=palette_values)
  if(!is.null(PlotTitle)){
    PosPlot<-PosPlot + ggtitle(PlotTitle) + theme(plot.title=element_markdown(size=8))
  }  
  if(!PlotLegend){
    PosPlot<-PosPlot+theme(legend.position="none")
  }
  return(PosPlot)
}

PlotCumCor<-function(PosDat,PredCol,palette_values,PlotPredicted=TRUE,
                     PlotLegend=TRUE,PlotTitle=NULL){
  # observed data
  PosDat$CumPres <- round(PosDat$CumPres) # round to only integers - decimals result from ambiguities in scoring  
  ObsTableCumPres<-PosDat %>% group_by(CumPres,pos) %>% summarise(Pres=mean(preserved),N=length(preserved))
  ObsTableCumPres$CumPres<-as.factor(ObsTableCumPres$CumPres)
  ObsTableCumPres<-ObsTableCumPres[ObsTableCumPres$N>10,]
  
  if(PlotPredicted){
    PredTableCumPres<-PosDat %>% group_by(CumPres,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(mean))
    PredTableN<-PosDat %>% group_by(CumPres,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(length))
    names(PredTableCumPres)[length(PredTableCumPres[1,])]<-"PredPres"
    names(PredTableN)[length(PredTableN[1,])]<-"N"
    PredTableCumPres$CumPres<-as.factor(PredTableCumPres$CumPres)
    PredTableCumPres<-PredTableCumPres[PredTableN$N>10,]
  }
  
  # plot observed values
  CorPlot<-ggplot(ObsTableCumPres,aes(x=pos,y=Pres,
                                      group=CumPres,
                                      shape=CumPres,
                                      color=CumPres))+
    geom_point(size=2)+geom_line(size=1,linetype="solid")
  
  # plot predicted values
  if(PlotPredicted){
    CorPlot<-CorPlot+geom_point(data=PredTableCumPres,
                                aes(x=pos,y=PredPres,
                                    group=CumPres,
                                    shape=CumPres,
                                    color=CumPres),
                                size=2)+
      geom_line(data=PredTableCumPres,aes(x=pos,y=PredPres,
                                          group=CumPres,
                                          color=CumPres),
                linetype="dashed",size=1.1)
  }
  
  CorPlot<-CorPlot + scale_shape_discrete(name="Previous correct")
  CorPlot<-CorPlot + scale_linetype_discrete(name="Previous correct")
  CorPlot<-CorPlot + scale_x_continuous(name="Word position",
                                        breaks=c(1,2,3,4,5,6,7,8,9,10))  +
    theme(axis.title.x = element_text(size=8))
  CorPlot<-CorPlot + scale_y_continuous(name="p(preserved)")  +
    theme(axis.title.y = element_text(size=8))
  CorPlot <- CorPlot + scale_color_manual(values=palette_values)
  if(!is.null(PlotTitle)){
    CorPlot<-CorPlot + ggtitle(PlotTitle) + + theme(plot.title=element_markdown(size=8))
  }
  if(!PlotLegend){
    CorPlot<-CorPlot+theme(legend.position="none")
  }
  return(CorPlot)
}

PlotCumErr<-function(PosDat,PredCol,palette_values,
                     PlotPredicted=TRUE,PlotLegend=TRUE,PlotTitle=NULL){
  # observed data
  PosDat$CumErr<-round(PosDat$CumErr) # so only integers are plotted (decimals are the result of ambiguities in scoring)
  ObsTableCumErr<-PosDat %>% group_by(CumErr,pos) %>% summarise(Pres=mean(preserved),N=length(preserved))
  ObsTableCumErr$CumErr<-as.factor(ObsTableCumErr$CumErr)
  ObsTableCumErr<-ObsTableCumErr[ObsTableCumErr$N>10,]
  
  # PredTableCumErr<-PosDat %>% group_by(CumErr,pos) %>% summarise(PredPres=mean(PosDat[,PredCol]),N=length(PosDat[,PredCol]))
  
  if(PlotPredicted){
    PredTableCumErr<-PosDat %>% group_by(CumErr,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(mean))
    PredTableN<-PosDat %>% group_by(CumErr,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(length))
    names(PredTableCumErr)[length(PredTableCumErr[1,])]<-"PredPres"
    names(PredTableN)[length(PredTableN[1,])]<-"N"
    PredTableCumErr$CumErr<-as.factor(PredTableCumErr$CumErr)
    PredTableCumErr<-PredTableCumErr[PredTableN$N>10,]
  }
  
  # plot observed values
  ErrPlot<-ggplot(ObsTableCumErr,aes(x=pos,y=Pres,
                                     group=CumErr,
                                     shape=CumErr,
                                     color=CumErr))+
    geom_point(size=2)+geom_line(size=1,linetype="solid")
  
  # plot predicted values
  if(PlotPredicted){
    ErrPlot<-ErrPlot+geom_point(data=PredTableCumErr,
                                aes(x=pos,y=PredPres,
                                    group=CumErr,
                                    shape=CumErr,
                                    color=CumErr),
                                size=2)+
      geom_line(data=PredTableCumErr,
                aes(x=pos,y=PredPres,
                    group=CumErr,
                    linetype=CumErr,
                    color=CumErr),linetype="dashed",
                size=1.1)
  }
  
  ErrPlot<-ErrPlot + scale_shape_discrete(name="Previous error")
  ErrPlot<-ErrPlot + scale_linetype_discrete(name="Previous error")
  ErrPlot<-ErrPlot + scale_x_continuous(name="Word position",
                                        breaks=c(1,2,3,4,5,6,7,8,9,10))  +
    theme(axis.title.x = element_text(size=8))
  ErrPlot<-ErrPlot + scale_y_continuous(name="p(preserved)")  +
    theme(axis.title.y = element_text(size=8))
  ErrPlot<-ErrPlot + scale_color_manual(values=palette_values)
  if(!is.null(PlotTitle)){
    ErrPlot<-ErrPlot + ggtitle(PlotTitle) + theme(plot.title=element_markdown(size=8))
  }  
  if(!PlotLegend){
    ErrPlot<-ErrPlot+theme(legend.position="none")
  }
  return(ErrPlot)
}

PlotFreq<-function(PosDat,PredCol,palette_values,
                   PlotPredicted=TRUE,PlotLegend=TRUE,PlotTitle=NULL){
  # observed data
  ObsTableFreq<-PosDat %>% group_by(freq_bin,pos) %>% summarise(Pres=mean(preserved),N=length(preserved))
  ObsTableFreq$freq_bin<-as.factor(ObsTableFreq$freq_bin)
  ObsTableFreq<-ObsTableFreq[ObsTableFreq$N>10,]
  
  # PredTableFreq<-PosDat %>% group_by(freq_bin,pos) %>% summarise(PredPres=mean(PosDat[,PredCol]),N=length(PosDat[,PredCol]))
  
  if(PlotPredicted){
    PredTableFreq<-PosDat %>% group_by(freq_bin,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(mean))
    PredTableN<-PosDat %>% group_by(freq_bin,pos) %>% summarise_at(.vars = colnames(.)[PredCol],.funs=c(length))
    names(PredTableFreq)[length(PredTableFreq[1,])]<-"PredPres"
    names(PredTableN)[length(PredTableN[1,])]<-"N"
    PredTableFreq$freq_bin<-as.factor(PredTableFreq$freq_bin)
    PredTableFreq<-PredTableFreq[PredTableN$N>10,]
  }
  
  # plot observed values
  FreqPlot<-ggplot(ObsTableFreq,aes(x=pos,y=Pres,
                                    group=freq_bin,
                                    shape=freq_bin,
                                    color=freq_bin))+
    geom_point(size=2)+geom_line(size=1,linetype="solid")
  
  # plot predicted values
  if(PlotPredicted){
    FreqPlot<-FreqPlot+geom_point(data=PredTableFreq,
                                  aes(x=pos,y=PredPres,
                                      group=freq_bin,
                                      shape=freq_bin,
                                      color=freq_bin),
                                  size=2)+
      geom_line(data=PredTableFreq,
                aes(x=pos,y=PredPres,
                    group=freq_bin,
                    color=freq_bin),linetype="dashed",size=1.1)
  }
  
  FreqPlot<-FreqPlot + scale_shape_discrete(name="Frequency bin")
  FreqPlot<-FreqPlot + scale_linetype_discrete(name="Frequency bin")
  FreqPlot<-FreqPlot + scale_x_continuous(name="Word position",
                                          breaks=c(1,2,3,4,5,6,7,8,9,10)) +
    theme(axis.title.x = element_text(size=8))
  FreqPlot<-FreqPlot + scale_y_continuous(name="p(preserved)") +
    theme(axis.title.y = element_text(size=8))
  FreqPlot<-FreqPlot + scale_color_manual(values=palette_values)
  if(!is.null(PlotTitle)){
    FreqPlot<-FreqPlot + ggtitle(PlotTitle)+ + theme(plot.title=element_markdown(size=8))
  }  
  if(!PlotLegend){
    FreqPlot<-FreqPlot+theme(legend.position="none")
  }
  return(FreqPlot)
}

ConvertDominanceResultToTable<-function(DA.Result){
  CA<-DA.Result$contribution.average
  CAunl<-unlist(CA)
  NumVars<-length(CA$r2.m)
  VarNames<-names(CA$r2.m)
  ResultMatrix<-matrix(CAunl,ncol=NumVars,byrow=TRUE)
  ResultDF<-as.data.frame(ResultMatrix)
  rownames(ResultDF)<-c("McFadden","SquaredCorrelation","Nagelkerke","Estrella")
  colnames(ResultDF)<-VarNames
  return(ResultDF)
}

ExtractDeviance<-function(Model){
  DevianceValue<-vector("numeric",length=1)
  DevianceValue[1]<-Model$deviance
  StringFormula<-as.character(Model$formula)
  names(DevianceValue)<-paste0(StringFormula[3:length(StringFormula)])
  return(DevianceValue)
}

GetDevianceSet<-function(ModelSet,MyData){
  ModelSetResults<-TestModels(ModelSet,MyData,PrintOutput=FALSE)
  DevianceSet<-sapply(ModelSetResults$ModelResult,ExtractDeviance)
  NullDeviance<-ModelSetResults$ModelResult[[1]]$null.deviance
  DevianceSetWithNull <- c(DevianceSet,NullDeviance)
  names(DevianceSetWithNull) <- c(names(DevianceSet),"null")
  OrderedDevianceSet<-DevianceSetWithNull[order(DevianceSetWithNull)]
  deviance_explained=OrderedDevianceSet["null"]-OrderedDevianceSet
  percent_explained=((OrderedDevianceSet["null"]-OrderedDevianceSet)/
                       OrderedDevianceSet["null"])*100
  percent_of_explained_deviance <- (deviance_explained / deviance_explained[1])*100
  increment_in_explained<-c(rev(diff(rev(percent_of_explained_deviance))),0)
  percent_of_explained_deviance["null"]<-NA
  DevianceDifferences<-c(0,diff(OrderedDevianceSet))
  DevianceSetDF<-data.frame(model=names(OrderedDevianceSet),
                            deviance=OrderedDevianceSet,
                            deviance_explained=deviance_explained,
                            percent_explained=percent_explained,
                            percent_of_explained_deviance=percent_of_explained_deviance,
                            increment_in_explained=increment_in_explained
  )
  return(DevianceSetDF)
}


# get observed and predicted values (from model formula model) for
# three ways of looking at/summarising the data: len x pos,
# cumulative error x pos and cumulative correct x pos. 
# the output is summarised data for observed values based on len x pos,
# cumulative error and cumulative correct and then output for model predicted
# values for len x pos, cumulative error and cumulative correct. 
# The format is two lists that have three lists each of p(preserved) values
# observed list -> len x pos list, cum error list, cum correct list
# predicted list > len x pos list, cum error list, cum correct list
# model formula (reported back for bookkeeping purposes)
# N_cutoff removes cells that have low N (noisy/poor estimates)
get_preserved_values<-function(model_formula, data_for_model, N_cutoff){
  obs_len_pos_table <- data_for_model %>% group_by(stimlen,pos) %>%
    summarise(preserved=mean(preserved),N=n()) %>% filter(N>N_cutoff)
  obs_cumerr_pos_table <- data_for_model %>% group_by(CumErr,pos) %>%
    summarise(preserved=mean(preserved),N=n()) %>% filter(N>N_cutoff)
  obs_cumcor_pos_table <- data_for_model %>% group_by(CumPres,pos) %>%
    summarise(preserved = mean(preserved),N=n()) %>% filter(N>N_cutoff)
  obs_preserved <- list(len_pos=obs_len_pos_table$preserved, 
                        cumerr=obs_cumerr_pos_table$preserved,
                        cumcor=obs_cumcor_pos_table$preserved)
  N_values<-list(len_pos_N=obs_len_pos_table$N,
                 cumerr_N=obs_cumcor_pos_table$N,
                 cumcor_N=obs_cumcor_pos_table$N)
  
  current_model <- glm(model_formula,data=data_for_model,family="binomial")
  model_deviance <- current_model$deviance
  data_for_model$fitted_values <- fitted(current_model)
  
  model_len_pos_table <- data_for_model %>% group_by(stimlen,pos) %>%
    summarise(fitted_values=mean(fitted_values),N=n()) %>% filter(N>N_cutoff)
  model_cumerr_pos_table <- data_for_model %>% group_by(CumErr,pos) %>%
    summarise(fitted_values=mean(fitted_values),N=n()) %>% filter(N>N_cutoff)
  model_cumcor_pos_table <- data_for_model %>% group_by(CumPres,pos) %>%
    summarise(fitted_values = mean(fitted_values),N=n()) %>% filter(N>N_cutoff)
  model_fitted <- list(len_pos=model_len_pos_table$fitted_values, 
                       cumerr=model_cumerr_pos_table$fitted_values,
                       cumcor=model_cumcor_pos_table$fitted_values)
  preserved_list<-list(obs_preserved, model_fitted, N_values, model_formula, model_deviance)
  return(preserved_list)
}

# This function takes one list of observed values (e.g. the first element of output
# of the function above), one list of model_fitted values (e.g. the second element
# of the output of the function above), values for the null model and the model
# formula for the model being tested (the third element of the function above)
# The function calculates the SS of the residuals from the model, the SS of the
# residuals from the null model and calculates the proportion of the null SS
# that is accounted for by the model
#
# IMPORTANT: This is based on an _unweighted_ set of summary values, and so results
# won't match exactly those from weighted data (e.g. based on deviance values).  This
# does, however, give a ballpark figure of what percentage of variance in the summarised data
# (e.g. data that might be plotted) that is accounted for by the model
# raw figures (as opposed to something based on a data summary) are sometimes hard
# to interpret because you don't expect to be able to predict trial-by-trial outcomes
# very accurately, but you do expect to do better on summarised 
# data (but then -- which summary?).  We use a summary that is defined by our theoretical
# concerns.
SSE_for_model_unweighted<-function(obs_preserved, model_fitted_values, null_fitted_values, model_formula){
  len_pos_residuals<-model_fitted_values$len_pos - obs_preserved$len_pos
  len_pos_resid_sqr<-(len_pos_residuals^2)
  cumerr_residuals<-model_fitted_values$cumerr - obs_preserved$cumerr
  cumerr_resid_sqr<-(cumerr_residuals^2)
  cumcor_residuals<-model_fitted_values$cumcor - obs_preserved$cumcor
  cumcor_resid_sqr<-(cumcor_residuals^2)
  residual_sse <- sum(len_pos_resid_sqr) +
    sum(cumerr_resid_sqr) +
    sum(cumcor_resid_sqr)  # what is left after current model predicts
  
  null_len_pos_resid<-null_fitted_values$len_pos - obs_preserved$len_pos
  null_len_pos_resid_sqr <- null_len_pos_resid^2
  null_cumerr_resid <- null_fitted_values$cumerr - obs_preserved$cumerr
  null_cumerr_resid_sqr <- null_cumerr_resid^2
  null_cumcor_resid <- null_fitted_values$cumcor - obs_preserved$cumcor
  null_cumcor_resid_sqr <- null_cumcor_resid^2
  all_sse <- sum(null_len_pos_resid_sqr) +
    sum(null_cumerr_resid_sqr) +
    sum(null_cumcor_resid_sqr)    
  
  # residuals from null model - residuals from current model = what is accounted for by the current model (model_ss)
  model_ss <- all_sse - residual_sse
  proportion_sse_acc_for<-model_ss/all_sse # % of all variation in the null model that is accounted for by the current model
  
  return_list<-list(all_sse=all_sse, model_ss=model_ss, residual_sse=residual_sse, proportion_accounted_for = proportion_sse_acc_for, model_formula=model_formula)
  return(return_list)
}


# weighted version
SSE_for_model_weighted<-function(obs_preserved, model_fitted_values, null_fitted_values, N_values, model_formula, model_deviance){
  len_pos_residuals<-model_fitted_values$len_pos - obs_preserved$len_pos
  len_pos_resid_sqr<-(len_pos_residuals^2)*(N_values$len_pos_N)
  cumerr_residuals<-model_fitted_values$cumerr - obs_preserved$cumerr
  cumerr_resid_sqr<-(cumerr_residuals^2)*(N_values$cumerr_N)
  cumcor_residuals<-model_fitted_values$cumcor - obs_preserved$cumcor
  cumcor_resid_sqr<-(cumcor_residuals^2)*(N_values$cumcor_N)
  residual_sse <- (sum(len_pos_resid_sqr) / sum(N_values$len_pos_N) )+
    (sum(cumerr_resid_sqr) / sum(N_values$cumerr_N)) +
    (sum(cumcor_resid_sqr) / sum(N_values$cumerr_N))  # what is left after current model predicts
  
  null_len_pos_resid<-null_fitted_values$len_pos - obs_preserved$len_pos
  null_len_pos_resid_sqr <- (null_len_pos_resid^2) *
    (N_values$len_pos_N)
  null_cumerr_resid <- null_fitted_values$cumerr - obs_preserved$cumerr
  null_cumerr_resid_sqr <- (null_cumerr_resid^2) *
    (N_values$len_pos_N)
  null_cumcor_resid <- null_fitted_values$cumcor - obs_preserved$cumcor
  null_cumcor_resid_sqr <- (null_cumcor_resid^2) *
    (N_values$len_pos_N)
  all_sse <- (sum(null_len_pos_resid_sqr) / sum(N_values$len_pos_N)) +
    (sum(null_cumerr_resid_sqr) / sum(N_values$cumerr_N)) +
    (sum(null_cumcor_resid_sqr) / sum(N_values$cumerr_N))   
  
  # residuals from null model - residuals from current model = what is accounted for by the current model (model_ss)
  model_ss <- all_sse - residual_sse
  proportion_sse_acc_for<-model_ss/all_sse # % of all variation in the null model that is accounted for by the current model
  
  return_list<-list(all_sse=all_sse, model_ss=model_ss, residual_sse=residual_sse, proportion_accounted_for = proportion_sse_acc_for, model_formula=model_formula, model_deviance = model_deviance)
  return(return_list)
}


# some other accessor or checking functions
get_proportion_accounted_for <- function(sse_result_list){
  return(sse_result_list$proportion_accounted_for)
}

get_model_formula <- function(sse_result_list){
  return(sse_result_list$model_formula)
}

get_deviance <- function(sse_result_list){
  return(sse_result_list$model_deviance)
}

# plot_table<-function(summary_table,line_var,dep_var){
#   summary_table<-summary_table[summary_table$N > 10,]
#   summary_table[line_var] <- factor(pull(summary_table,line_var))
#   summary_table["pos"]<- factor(pull(summary_table,pos))
#   
#   table_plot <- ggplot(summary_table,aes(x=pos,
#                                          y=.data[[dep_var]],
#                                          group=.data[[line_var]],
#                                          color=.data[[line_var]],
#                                          shape=.data[[line_var]],
#                                          linetype = .data[[line_var]]),color="black") + 
#     geom_point(size=3) +
#     geom_line()
#   return(table_plot)
# }


compare_SS_accounted_for<-function(ModelList,NullModel,MyData,N_cutoff){
  model_preserved_values<-list()
  num_models<-length(ModelList)
  
  # null model
  null_fitted_values <- get_preserved_values(NullModel,MyData,N_cutoff)
  
  sse_results_list <- list()
  
  for(i in seq(1,num_models)){
    model_preserved_values[[i]]<-get_preserved_values(ModelList[i],MyData,N_cutoff)
    current_pres_values <- model_preserved_values[[i]]
    sse_results_list[[i]]<-SSE_for_model_weighted(current_pres_values[[1]], # observed values
                                                  current_pres_values[[2]], # model predicted values
                                                  null_fitted_values[[2]],  # null model predicted values
                                                  current_pres_values[[3]], # N values
                                                  current_pres_values[[4]], # model formula
                                                  current_pres_values[[5]]) # model deviance (for checking)
  }
  return(sse_results_list)
}

difference_matrix<-function(value_vector){
  num_values<-length(value_vector)
  difference_matrix <- matrix(nrow=num_values, ncol=num_values)
  for(i in seq(1,num_values-1)){ # last row is 0 diagonal only
    difference_matrix[i,i]<-0
    for(j in seq(i+1,num_values)){ # start at i+1 - upper triangle/lower triangle 
      difference_matrix[j,i] <- value_vector[j] - value_vector[i]
      difference_matrix[i,j] <- value_vector[i] - value_vector[j]
    }
  }
  difference_matrix[num_values,num_values]<-0
  return(difference_matrix)
}

sse_results_table<-function(list_of_sse_results){
  num_models <- length(list_of_sse_results)
  proportion_accounted_for <- sapply(list_of_sse_results,get_proportion_accounted_for) 
  sorted_order <- order(proportion_accounted_for)
  sorted_results <- list_of_sse_results[sorted_order]  
  proportion_accounted_for_sorted <- sapply(sorted_results,get_proportion_accounted_for) 
  model_deviance_sorted <- sapply(sorted_results,get_deviance)
  differences <- difference_matrix(proportion_accounted_for_sorted)
  model_formulas_sorted <- sapply(sorted_results,get_model_formula)
  results_table<-data.frame(model=model_formulas_sorted)
  results_table$p_accounted_for <- proportion_accounted_for_sorted
  results_table$model_deviance <- model_deviance_sorted
  for(i in seq(1,num_models)){ # columns
    col_name <- sub("preserved ~ ","diff_",results_table$model[i])
    results_table[col_name]<-differences[,i]
  }
  return(results_table)
}


