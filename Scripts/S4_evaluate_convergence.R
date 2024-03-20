##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 4: Evaluates convergence of the fitted models
## Author: Maria Rakka
##----------------------------------------------------------------------------


## Libraries
library(Hmsc)
library(colorspace)
library(vioplot)

## Set directories
setwd("...")
localDir = "."
modelDir = file.path(localDir, "2_Models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

##Make the script reproducible
set.seed(1)

## Set parameters for loop

showBeta = TRUE  # Show convergence for beta-parameters
showGamma = TRUE # Show convergence for gamma-parameters
showOmega = TRUE # Show convergence for omega-parameters
maxOmega = 100   # Show convergence for 100 randomly selected species pairs
showRho = TRUE   # Show convergence for rho-parameters
showAlpha = TRUE # Show convergence for alpha-parameters

if(is.null(showBeta)) showBeta = TRUE
if(is.null(showGamma)) showGamma = FALSE
if(is.null(showOmega)) showOmega = FALSE
if(is.null(maxOmega)) maxOmega = 50
if(is.null(showRho)) showRho = FALSE
if(is.null(showAlpha)) showAlpha = FALSE

## Define how many samples, thinning and chains you used so that the code
## picks the right fitted models

samples_list = c(5,250)
thin_list = c(1,100)
nst = length(thin_list)
nChains = 4

## Produce graphs for convergence in pdf file

text.file = file.path(resultDir,"/MCMC_convergence.txt")
cat("MCMC Convergennce statistics\n\n",file=text.file,sep="")

ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL  
ma.rho = NULL
na.rho = NULL
Lst = 1
while(Lst <= nst){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){
    load(filename)
    cat(c("\n",filename,"\n\n"),file=text.file,sep="",append=TRUE)
    nm = length(models)
    for(j in 1:nm){
      mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      nr = models[[j]]$nr
      cat(c("\n",names(models)[j],"\n\n"),file=text.file,sep="",append=TRUE)
      if(showBeta){
        psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\nbeta\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.beta)){
          ma.beta = psrf[,1]
          na.beta = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.beta = cbind(ma.beta,psrf[,1])
          if(j==1){
            na.beta = c(na.beta,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.beta = c(na.beta,"")
          }
        }
      }
      if(showGamma){
        psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\ngamma\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.gamma)){
          ma.gamma = psrf[,1]
          na.gamma = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.gamma = cbind(ma.gamma,psrf[,1])
          if(j==1){
            na.gamma = c(na.gamma,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.gamma = c(na.gamma,"")
          }
        }
      }
      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat("\nrho\n\n",file=text.file,sep="",append=TRUE)
        cat(psrf[1],file=text.file,sep="\n",append=TRUE)
      }
      if(showOmega & nr>0){
        cat("\nomega\n\n",file=text.file,sep="",append=TRUE)
        for(k in 1:nr){
          cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
          tmp = mpost$Omega[[k]]
          z = dim(tmp[[1]])[2]
          if(z > maxOmega){
            sel = sample(1:z, size = maxOmega)
            for(i in 1:length(tmp)){
              tmp[[i]] = tmp[[i]][,sel]
            }
          }
          psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
          tmp = summary(psrf)
          cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
          if(is.null(ma.omega)){
            ma.omega = psrf[,1]
            na.omega = paste0(as.character(thin),",",as.character(samples))
          } else {
            ma.omega = cbind(ma.omega,psrf[,1])
            if(j==1){
              na.omega = c(na.omega,paste0(as.character(thin),",",as.character(samples)))
            } else {
              na.omega = c(na.omega,"")
            }
          }
        }
      }
      if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models[[j]]$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n",file=text.file,sep="\n",append=TRUE)
            cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=text.file,sep="\n",append=TRUE)            
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

pdf(file= file.path(resultDir,"/MCMC_convergence.pdf"))
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}
dev.off()
