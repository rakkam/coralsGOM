##--------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 1: Creates a 3D trait space including all the coral traits
## Author: Maria Rakka
##--------------------------------------------------------------------------


## Libraries
library(tidyverse)
library(mFD)

## Make script reproducible
set.seed(1)

## Set directories
setwd("...")
localDir = "."
rawdataDir = file.path(localDir, "1_RawData")
modelDir = file.path(localDir, "2_Models")

## Read model data
load(file.path(rawdataDir,"4model_data_1km.RData"))

trait_dat<-trait_dat %>% 
  select(Aragonite,Calcite,Scleroprotein,Polyp_diameter,Height) %>% 
  as.data.frame()

rownames(trait_dat)<-names(model_cordat)

## Create a 3D trait space

traits_cat <- data.frame(names(trait_dat), 
                         c("F","F","F","Q","Q"),
                         c("skeltrait","skeltrait","skeltrait",NA,NA))
colnames(traits_cat) <- c("trait_name", "trait_type","fuzzy_name")
traits_cat

dist_gen <- mFD::funct.dist(
  sp_tr         = trait_dat,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

summary(as.matrix(dist_gen))

round(dist_gen,3)

## Evaluate the quality of the constructed space

quality_fspaces_fruits <- mFD::quality.fspaces(
  sp_dist             = dist_gen,
  fdendro             = "average",
  maxdim_pcoa         = 9,
  deviation_weighting = c("absolute", "squared"),
  fdist_scaling       = c(TRUE, FALSE))

quality_fspaces_fruits$"details_fspaces"$"dendro" %>%
  as.dendrogram() %>%
  dendextend::plot_horiz.dendrogram(side = TRUE)

round(quality_fspaces_fruits$"quality_fspaces", 3)

quality.fspaces.plot(
  fspaces_quality            = quality_fspaces_fruits,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

sp_faxes_coord_fruits <- quality_fspaces_fruits$"details_fspaces"$"sp_pc_coord"

fruits_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = trait_dat, 
  sp_faxes_coord = sp_faxes_coord_fruits[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

fruits_tr_faxes$"tr_faxes_plot"+
  geom_text(aes(x=0.1,y=50),"cor")

sp_faxes_coord_fruits <- quality_fspaces_fruits$"details_fspaces"$"sp_pc_coord"

trait_axes<-as.data.frame(sp_faxes_coord_fruits) %>%
  select(PC1,PC2,PC3)

trait_axes_writable<-trait_axes %>% 
  rownames_to_column("Genus")

## Save trait table
write_csv(trait_axes,file=file.path(rawdataDir,"AxesTraits_1km.csv"))

fruits_tr_faxes$"tr_faxes_stat"[which(fruits_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

## Create final graphs
sig.axis<-fruits_tr_faxes$tr_faxes_stat %>% 
  mutate(sig=ifelse(p.value<0.05,"yes","no")) %>% 
  select(trait,axis,value,sig)

sig.labs<-sig.axis %>% 
  filter(sig=="yes") %>% 
  mutate(pcoval=case_when(axis=="PC1"~0.45,
                          axis=="PC2"~0.35,
                          TRUE~0.2),
         traitval=case_when(trait=="Height" & axis=="PC3"~6,
                            trait=="Polyp_diameter"& axis=="PC3"~3,
                            trait=="Polyp_diameter" & axis %in%c("PC1","PC2")~0.5,
                            trait=="Height"& axis %in%c("PC1","PC2")~3,
                            TRUE~25),
         value=as.character(paste("~R^{2}==",round(value,2),sep="")))%>% 
  mutate(trait=fct_relevel(fct_recode(trait,'Polyp diameter'="Polyp_diameter"),
                           "Aragonite","Calcite","Scleroprotein"))

cbind(trait_dat,trait_axes) %>%
  pivot_longer(cols=c("Aragonite","Calcite","Scleroprotein","Polyp_diameter","Height"),
               names_to = "trait",values_to = "traitval") %>% 
  pivot_longer(cols=c("PC1","PC2","PC3"),names_to="axis",values_to = "pcoval") %>%
  left_join(sig.axis,by=c("trait","axis")) %>% 
  mutate(trait=fct_relevel(fct_recode(trait,'Polyp diameter'="Polyp_diameter"),
                           "Aragonite","Calcite","Scleroprotein")) %>%
  ggplot(aes(x=traitval,y=pcoval))+
  geom_point(aes(col=sig),size=2)+
  geom_text(data=sig.labs,
            mapping = aes( y = +Inf, label = value),
            vjust   = 1.8,colour="grey20",size=3.5,parse=TRUE)+
  facet_grid(vars(axis),vars(trait),scales="free",switch="both")+
  theme_light()+
  scale_colour_manual(values=c("grey","#6665DD"))+
  theme(axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour="black",face="bold"),
        legend.position="none",
        strip.placement = "outside")

