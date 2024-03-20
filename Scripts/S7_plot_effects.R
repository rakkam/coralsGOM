##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 7: Plot model effects 
## Author: Maria Rakka
##----------------------------------------------------------------------------

#libraries
library(tidyverse)
library(Hmsc)

## Prepare directories
setwd("...")
localDir = "."
modelDir = file.path(localDir, "2_Models")
outputDir=file.path(localDir,"3_Outputs")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

## Prepare directories
load(file.path(modelDir,"models_thin_100_samples_250_chains_4.RData"))

## Calculate total and marginal effects
m<-models$pres.abs.quad
covariates <- all.vars(m$XFormula)
cov.names<-c("Salinity","Temperature","Current velocity",
             "TPI","Mud (%)","Aspect","Sampling effort")
myspecies_names<-paste(models$pres.abs.quad$spNames,"sp",sep=" ")

#Make functions to calculate effects for each variable
calc_total_eff<-function(myvariable){
  Gradient = constructGradient(m,focalVariable = myvariable)
  predY = predict(m, Gradient=Gradient, expected = TRUE) 
return(list(mygrads=Gradient,mypreds=predY))}

calc_marg_eff<-function(myvariable){
  Gradient = constructGradient(m,focalVariable = myvariable,non.focalVariables = 1)
  predY = predict(m, Gradient=Gradient, expected = TRUE) 
  return(list(mygrads=Gradient,mypreds=predY))}

#Apply function to all variables
preps<-lapply(covariates,calc_total_eff)
preps_marg<-lapply(covariates,calc_marg_eff)

# Save data to use in the future without estimating again

save(preps,preps_marg,file=file.path(outputDir, "curv_dat.RData"),sep="/")

##Make graphs

#Generate graphs for environmental variables
pdf(file.path(resultDir,"1.Total effects genera.pdf"))
par(mfrow=c(4,4),mar=c(2.5,2,2,1),mgp=c(2,0.6,0),oma = c(2, 0, 6, 0))
for(i in 1:length(covariates)){
  for (j in 1:30){
    plotGradient(m, preps[[i]]$mygrads, pred=preps[[i]]$mypreds,
                 measure="Y",index=j, showData = TRUE,
                 yshow = 0,
                 cicol=rgb(0.63,0.768,0.737,alpha=.3),
                 showPosteriorSupport=FALSE,main=myspecies_names[[j]],
                 xlabel="",cex.axis=0.7)
    mtext(side=1, text=cov.names[[i]], line=1.75,cex= 0.7)
}
}

dev.off()

pdf(file.path(resultDir,"2.Marginal effects genera.pdf"))
par(mfrow=c(4,4),mar=c(2.5,2,2,1),mgp=c(2,0.6,0),oma = c(2, 0, 6, 0))
for(i in 1:length(covariates)){
  for (j in 1:30){
    plotGradient(m, preps_marg[[i]]$mygrads, pred=preps_marg[[i]]$mypreds,
                 measure="Y",index=j, showData = TRUE,
                 yshow = 0,
                 cicol=rgb(0.63,0.768,0.737,alpha=.3),
                 showPosteriorSupport=FALSE,main=myspecies_names[[j]],
                 xlabel="",cex.axis=0.7)
    mtext(side=1, text=cov.names[[i]], line=1.75,cex= 0.7)
  }
}
dev.off()


#Generate graphs for traits
mytraits<-c("intercept",all.vars(m$TrFormula))

pdf(file.path(resultDir,"3.Total effects traits.pdf"))
par(mfrow=c(4,3),mar=c(2.5,2,2,1),mgp=c(2,0.6,0),oma = c(2, 0, 6, 0))
for(i in 1:length(covariates)){
  for (j in 2:length(mytraits)){
    plotGradient(m, preps[[i]]$mygrads, pred=preps[[i]]$mypreds,
                 measure="T",index=j, showData = TRUE,
                 cicol=rgb(0.63,0.768,0.737,alpha=.3),
                 showPosteriorSupport=FALSE,main=mytraits[j],
                 xlabel="",cex.axis=0.7)
    mtext(side=1, text=cov.names[[i]], line=1.75,cex= 0.7)
  }
}
dev.off()

pdf(file.path(resultDir,"4.Marginal effects traits.pdf"))
par(mfrow=c(4,3),mar=c(2.5,2,2,1),mgp=c(2,0.6,0),oma = c(2, 0, 6, 0))
for(i in 1:length(covariates)){
  for (j in 2:length(mytraits)){
    plotGradient(m, preps_marg[[i]]$mygrads, pred=preps_marg[[i]]$mypreds,
                 measure="T",index=j, showData = TRUE,
                 cicol=rgb(0.63,0.768,0.737,alpha=.3),
                 showPosteriorSupport=FALSE,main=mytraits[j],
                 xlabel="",cex.axis=0.7)
    mtext(side=1, text=cov.names[[i]], line=1.75,cex= 0.7)
  }
}
dev.off()
