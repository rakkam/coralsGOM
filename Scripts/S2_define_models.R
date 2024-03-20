##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 2: Takes original datafiles and defines the model structure
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Libraries
library(tidyverse)
library(Hmsc)

## Make script reproducible
set.seed(1)

## Set directories
setwd("...")
localDir = "."
rawdataDir = file.path(localDir, "1_RawData")
modelDir = file.path(localDir, "2_Models")

## Read and prepare data

load(paste(rawdataDir,"4model_data_1km.RData",sep="/"))
final_model_traits<-read.csv(file.path(dataDir,"AxesTraits_1km.csv"))

# Prepare xy matrix for spatial random factors where row names are the IDs of 
# each quadrat, and the columns are the xy-coordinates

xycoords<-model_expdes %>% 
  dplyr::select(ID,x=Lon,y=Lat) 

xydat<-xycoords%>% 
  column_to_rownames("ID") %>% 
  as.matrix()

# Prepare the studyDesign as a dataframe 
exdes<-model_expdes%>%
  dplyr::select(ID,Dive)%>% 
  mutate_all(as_factor) %>% 
  as.data.frame()

envdat<-model_envdat %>% 
  dplyr::select(-c(SurveyID,Dive)) %>% 
  column_to_rownames("ID")

rownames(final_model_traits)<-names(model_cordat)

# Prepare coral data (Ydata) with presence/absence only, arrange in matrix
model_cordat_pa<-model_cordat %>% 
  mutate_all(~ifelse(.>0,1,0))

## Set up the model

# Define the environmental model through XFormula
XFormula=~sal_mean+poly(temp_mean,degree=2,raw=TRUE)+vel_mean+tpi20km+mud+aspect2000+sam_ef


# Define the trait model through TrFormula
TrFormula = ~ PC1+PC2+PC3

# Set up the random effects on the Dive level and quadrat level (spatial)
rL.loc=HmscRandomLevel(units = levels(exdes$Dive))
rL.quad = HmscRandomLevel(sData = xydat)

# Use the Hmsc model constructor to define a model
# We used probit distribution as appropriate for presence-absence data

m=Hmsc(Y=model_cordat_pa,
       XData=envdat,XFormula=XFormula,
       TrData = final_model_traits, TrFormula = TrFormula,
       distr="probit",
       studyDesign = exdes,
       ranLevels=list(Dive=rL.loc,"ID"=rL.quad))

# Save models in list (not necessary if you only have 1 model)

models=list(m)
names(models) = c("pres.abs.quad")
save(models, file = file.path(modelDir, "unfitted_models.RData"))

## Test if models fit without errors

for(i in 1:length(models)){
  print(i)
  sampleMcmc(models[[i]],samples=2)
}
 