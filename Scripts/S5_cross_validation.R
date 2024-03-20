##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 5: Cross-validation
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Please note that running this script can take more than 48h so it is
## preferable to run it on a cluster computer


## Libraries
library(Hmsc)

## Set the base directory 
setwd("...")
localDir = "."
modelDir = file.path(localDir, "2_Models")

## Make script reproducible
set.seed(1)

## Set number of folds and parallel computing

nfolds = 10 #change the number of CV-folds
nParallel = 4 #Set to 1 to disable parallel computing


if(is.null(nfolds)) nfolds = 2

samples_list = c(5,250)
thin_list = c(1,100)
nChains = 4

if(is.null(nParallel)) nParallel = nChains
Lst = 1
while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename.in = file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
  filename.out = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                          "_samples_", as.character(samples),
                                          "_chains_",as.character(nChains),
                                          "_nfolds_", as.character(nfolds),
                                          ".Rdata",sep = ""))
  if(file.exists(filename.out)){
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    print("model fit had been computed already")
  } else {
    if(file.exists(filename.in)){
      print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
      print(date())
      load(file = filename.in) #models
      nm = length(models)
      
      MF = list()
      MFCV = list()
      WAIC = list()
      
      for(mi in 1:nm){
        print(paste0("model = ",names(models)[mi]))
        m = models[[mi]]
        preds = computePredictedValues(m)
        MF[[mi]] = evaluateModelFit(hM=m, predY=preds)
        partition = createPartition(m, nfolds = nfolds) #USE column = ...
        preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
        MFCV[[mi]] = evaluateModelFit(hM=m, predY=preds)
        WAIC[[mi]] = computeWAIC(m)
      }
      names(MF)=names(models)
      names(MFCV)=names(models)
      names(WAIC)=names(models)
      
      save(MF,MFCV,WAIC,file = filename.out)
    }
  }
  Lst = Lst + 1
}

