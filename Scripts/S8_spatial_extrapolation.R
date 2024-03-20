##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 7: Evaluates the quantity of analogous conditions
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Libraries
library(tidyverse)
library(terra)
library(raster)
library(dsmextra)
library(magrittr)

## Set and prepare directories
setwd("...")
localDir = "."
rawdataDir = file.path(localDir, "1_RawData")
outputDir=file.path(localDir,"3_Outputs")

## Get data of environmental variables for the whole study area
load(file=file.path(rawdataDir,"misc/Data_for_spatial_extrapolation.RData"))

testdat<-all_dat %>% 
  rename(x=Lon,y=Lat)

## Get grid of study area to copy crs
env_dat<-raster(file.path(rawdataDir,"misc/depth.grd"))

## Get environmental data that were used to train the model

load(file.path(rawdataDir,"4model_data_1km.RData"))

traindat<-cbind(model_envdat,model_expdes) %>% 
  dplyr::select(sal_mean,temp_mean,vel_mean,tpi20km,aspect2000,mud,x=Lon,y=Lat)

## Perform extrapolation

#Set tibble options
options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)


# Set knitr options
knitr::opts_chunk$set(echo = TRUE)

# Define projected coordinate system
aftt_crs <- crs(env_dat)

# Define environmental covariates of interest
covariates.spermwhale <- c("sal_mean","temp_mean","vel_mean","tpi20km","aspect2000","mud")

testextrapolation <- compute_extrapolation(samples = traindat,
                                                  covariate.names = covariates.spermwhale,
                                                  prediction.grid = testdat,
                                                  coordinate.system = aftt_crs)
summary(testextrapolation)

#Visualize

extrap_results<-map_extrapolation(map.type = "extrapolation",
                  extrapolation.object = testextrapolation)

extrap_vars<-map_extrapolation(map.type = "mic",
                                  extrapolation.object = testextrapolation)


# Calculate Gower's distances and %N

mynearby <- compute_nearby(samples = traindat,
                                      prediction.grid = testdat,
                                      coordinate.system = aftt_crs,
                                      covariate.names = covariates.spermwhale,
                                      nearby = 1)

map_extrapolation(map.type = "nearby",
                  extrapolation.object = mynearby)


analogue_areas<-testextrapolation$data$all %>%
  mutate(myrows=1:n()) %>%
  filter(ExDet>0 & ExDet<1)
  #filter(mic==0)

analogue_dat<-testdat[analogue_areas$myrows,]

# Save analogous areas to use for predictions

save(analogue_dat,file=file.path(outputDir,"Analogue_areas_spatial_extrapolation.RData"))

## Visualize extrapolation results (Make raster map)

dataformap<-testextrapolation$data$all %>% 
   mutate(mycat=case_when(ExDet<0~"Univariate",
                        ExDet>1 ~"Combinatorial",
                        ExDet>0 & ExDet<1 ~ "Analogous",
                        TRUE~NA)) %>% 
  dplyr::select(x,y,mycat)

ggplot(dataformap)+
  geom_raster(aes(x=x,y=y,fill=mycat))+
  coord_equal()+
  theme_minimal()+
  xlab("Lon")+ylab("Lat")+
  theme(legend.title=element_blank())
