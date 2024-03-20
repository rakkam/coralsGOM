##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 8: Produces predictions of species occurrences over study area
## (only using analogous areas)
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Libraries
library(tidyverse)
library(Hmsc)
library(foreach)
library(doParallel)
library(terra)
library(ggpubr)

## Set and prepare directories
setwd("...")
localDir = "."
rawdataDir = file.path(localDir, "1_RawData")
modelDir = file.path(localDir, "2_Models")
outputDir = file.path(localDir, "3_Outputs")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

## Get analogue areas

load(paste(outputDir,"Analogue_areas_spatial_extrapolation.RData",sep="/"))

analogue_dat<-analogue_dat %>% 
  rename(Lon=x,Lat=y)
split_vectors <- split(analogue_dat, cut(seq(nrow(analogue_dat)), 7))

## Load fitted model and make predictions

load(paste(modelDir,"models_thin_100_samples_250_chains_4.RData",sep="/"))

# Make function to generate predictions and metrics

get_preds<-function(mydat){
  trial_env1<-mydat %>% 
    dplyr::select(temp_mean,sal_mean,vel_mean,aspect2000,mud, tpi20km,sam_ef,ID=cell_nb) %>% 
    column_to_rownames("ID") %>% 
    as.data.frame()
  
  xydat<-mydat %>% 
    dplyr::select(Lon,Lat) %>% 
    #column_to_rownames("ID") %>% 
    as.matrix()
  
  
  gradient_trial=prepareGradient(models$pres.abs.quad,
                                 XDataNew=trial_env1,
                                 sDataNew=list("ID"=xydat))
  predict_trial=predict(models$pres.abs.quad,Gradient=gradient_trial,predictEtaMean = TRUE)
  EpredY=Reduce("+",predict_trial)/length(predict_trial)
  getS = function(p){return(rowSums(p))}
  aS = simplify2array(lapply(X = predict_trial, FUN = getS))
  dim(aS)
  Srich = apply(aS, 1, mean)
  sdS = sqrt(apply(aS, 1, var))
  S5 = apply(aS > 5, 1, mean)
  S10 = apply(aS > 5, 1, mean)
  CWM = (EpredY %*% models$pres.abs.quad$Tr)/matrix(rep(Srich, models$pres.abs.quad$nt), ncol = models$pres.abs.quad$nt)
  
  return(data.frame(Srich,sdS,S5,S10,CWM))
}

# Perform function in parallel

cl <- parallel::makeCluster(7)
doParallel::registerDoParallel(cl)
rf <- foreach(mydat=split_vectors, .packages=c('Hmsc',"tidyverse")) %dopar%
  get_preds(mydat=mydat)
parallel::stopCluster(cl)

spatial_preds<-cbind(do.call(rbind,rf),
      Lon=analogue_dat$Lon,Lat=analogue_dat$Lat)

# Save predictions

save(spatial_preds,file=paste(outputDir,"Spatial_predictions.RData",sep="/"))

## Plot predictions
# Get a random raster to copy crs

env_layer<-raster(file.path(rawdataDir,"misc/depth.grd"))

sp_ric_dat<-spatial_preds %>% 
  dplyr::select(x=Lon,y=Lat,Srich,sdS,S5,PC1,PC2,PC3)

sp_ric_ras <- rast(sp_ric_dat, type="xyz")
crs(sp_ric_ras) <- as.character(crs(env_layer))


##Prepare basemap
library(ggOceanMapsData)
library(ggOceanMaps)
library(ggnewscale)

load(file=paste(rawdataDir,"misc/shapefiles_for_basemap.RData",=sep="/"))

##Maps for manuscript
#const_colors = rev(c( "#bf5b1d", "#e48043", "#d6cdcd", "#84a7c4", "#547a99"))
const_colors = rev(c( "#bf5b1d", "#e48043", "#fee391", "#7fcdbb", "#547a99"))
sric<-basemap(shapefiles=list(land = bs_land, glacier = NULL,bathy=bs_bathy),
              c(-73,-63, 38.5,45),bathymetry = TRUE,grid.col = NA,
              bathy.style = "poly_greys",land.col = "#eeeac4",bathy.alpha = 0.2,
              land.border.col = NA)+
  annotation_scale(location = "br",style="ticks",text_col="white",line_col="white") + 
  annotation_north_arrow(location = "tr", which_north = "true",height = unit(1, "cm"),
                         width = unit(1, "cm"),) +
  #  geom_point(data = model_expdes, aes(x = Lon, y = Lat), 
  #             color = "black",shape=21,size=2.5,fill="goldenrod2")+
  labs(x = NULL, y = NULL)+
  guides(fill = "none")+
  ggnewscale::new_scale_fill()+
  geom_raster(data=sp_ric_dat,aes(
    x = x, 
    y = y, 
    fill = Srich))+
  scale_fill_gradientn(colours=const_colors)+
  theme(legend.position="bottom")

sric_sd<-basemap(shapefiles=list(land = bs_land, glacier = NULL,bathy=bs_bathy),
                 c(-73,-63, 38.5,45),bathymetry = TRUE,grid.col = NA,
                 bathy.style = "poly_greys",land.col = "#eeeac4",bathy.alpha = 0.2,
                 land.border.col = NA)+
  annotation_scale(location = "br",style="ticks",text_col="white",line_col="white") + 
  annotation_north_arrow(location = "tr", which_north = "true",height = unit(1, "cm"),
                         width = unit(1, "cm"),) +
  #  geom_point(data = model_expdes, aes(x = Lon, y = Lat), 
  #             color = "black",shape=21,size=2.5,fill="goldenrod2")+
  labs(x = NULL, y = NULL)+
  guides(fill = "none")+
  ggnewscale::new_scale_fill()+
  geom_raster(data=sp_ric_dat,aes(
    x = x, 
    y = y, 
    fill = sdS))+
  scale_fill_gradientn(colours=const_colors)+
  theme(legend.position="bottom")

ggarrange(sric,sric_sd,nrow=1,labels=c("A","B"))