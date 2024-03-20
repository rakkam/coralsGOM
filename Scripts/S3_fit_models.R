##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 3: Takes unfit models and fits them by using MCMC chains
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Please note that running this script can take more than 48h so it is
## preferrable to run it on a cluster computer

## Libraries
library(Hmsc)

# Set the base directory 
setwd("...")
localDir = "."
modelDir = file.path(localDir, "2_Models")

nParallel = NULL #Default: nParallel = nChains

## Load unfitted models that were created with S2_define_models.R script

load(file=file.path(modelDir,"unfitted_models.RData"))

## Fit models with 4 chains, 250 samples and thinning =100
## Models with less samples and thinning were also run, to be able
## to inspect them in cases of errors, or if the cluster computer would abort
## the submitted job before completing it

nm = length(models)
samples_list = c(5,250)
thin_list = c(1,100)
nChains = 4

if(is.null(nParallel)) nParallel = nChains
Lst = 1
while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                             "_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             ".Rdata",sep = ""))
  if(file.exists(filename)){
    print("model had been fitted already")
  } else {
    print(date())
    for (mi in 1:nm) {
      print(paste0("model = ",names(models)[mi]))
      m = models[[mi]]
      m = sampleMcmc(m, samples = samples, thin=thin,
                     adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains,
                     nParallel = nParallel) 
      models[[mi]] = m
    }
    save(models,file=filename)
  }
  Lst = Lst + 1
}
