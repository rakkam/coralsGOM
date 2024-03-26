##--------------------------------------------------------------------------
## Gulf of Maine models
## Regions of common profile, graphs and maps
## First created 27 Jun 2023
## Last modified 27 Jun 2023
##--------------------------------------------------------------------------

library(tidyverse)
library("ggsci")
library(raster)
library(sf)
library(ggpubr)
library(Hmsc)

## Set directories
setwd("...")
localDir = "."
rawdataDir = file.path(localDir, "1_RawData")
modelDir = file.path(localDir, "2_Models")
resultDir = file.path(localDir, "results")
outputDir=file.path(localDir,"3_Outputs")
if (!dir.exists(resultDir)) dir.create(resultDir)

## Load metrics
load(file=paste("RCP_species_info.RData",sep="/"))

mypalette<-c("#ABC4AB",
                      "#855A5C",
                      "#6E92A1",
                      "#080E3A",
                      "#E2C044",
                      "#A87438",
                      "#ABB557",
                      "#27556C",
                      "#8A8E91")
                      
cwm_graph<-mymetrics %>% 
  mutate( zone=as.factor(case_when(clust==1~"A",
                                   clust==2~"B",
                                   clust==3~"C",
                                   clust==4~"D",
                                   clust==7~"E",
                                   clust==5~ "F",
                                   clust==6~ "G",
                                   clust==9~"H",
                                   clust==8~"I"))) %>% 
  mutate(zone=fct_relevel(zone,c("A","B","C","D","E","F","G","H","I"))) %>% 
  dplyr::select(PC1,PC2,PC3,zone) %>% 
  pivot_longer(cols=!zone,values_to = "CWM",names_to="metric") %>% 
  ggplot(aes(x=zone,y=CWM)) +
  geom_jitter(aes(color=zone), alpha=0.2)+
  geom_boxplot(colour="gray30",alpha=0.4,outlier.shape = NA)+
  facet_wrap(~metric,ncol=1)+
  scale_colour_manual(values=mypalette)+
  #scale_fill_npg()+
  theme_light()+
  theme(strip.background = element_rect(fill="white",colour="gray"),
        strip.text = element_text(colour="black"),
        legend.title=element_blank(),
        legend.key.size = unit(7,"points"),
        legend.position="top")+
  xlab("Assemblage")+
  guides(colour = guide_legend(nrow = 1, byrow = TRUE,
         override.aes = list(alpha = 1)))


#Create dummy data to get better legend (legend has low alpha and points appear
#Transparent

dummyData <- data.frame(value=rep(NA,36),
                        zone=rep(LETTERS[1:9],4),
                        metric=rep(c("Srich","fric","fdis","feve"),9))


paper_metrics<-mymetrics %>% 
  mutate( zone=as.factor(case_when(clust==1~"A",
                                   clust==2~"B",
                                   clust==3~"C",
                                   clust==4~"D",
                                   clust==7~"E",
                                   clust==5~ "F",
                                   clust==6~ "G",
                                   clust==9~"H",
                                   clust==8~"I"))) %>% 
  mutate(zone=fct_relevel(zone,c("A","B","C","D","E","F","G","H","I")))%>% 
  dplyr::select(Srich,fric,feve,fdis,fdiv,zone) %>% 
  pivot_longer(cols=!zone,values_to = "value",names_to="metric") %>% 
  mutate(metric=fct_relevel(metric,c("Srich","fric","fdis","feve"))) %>%
  filter(!metric=="fdiv") %>%
  mutate(metric=recode_factor(metric,Srich="Species richness",
                              fric="Functional richness",
                              fdis="Functional dispersion",
                              feve="Functional evenness")) %>%
  ggplot(aes(x=zone,y=value)) +
  geom_jitter(aes(color=zone), alpha=0.2)+
  geom_boxplot(colour="gray30",alpha=0.4,outlier.shape = NA)+
  facet_wrap(~metric,ncol=2,scales="free")+
  scale_colour_manual(values=mypalette)+
  theme_light()+
  theme(strip.background = element_rect(fill="white",colour="gray"),
        strip.text = element_text(colour="black"),
        legend.title=element_blank(),
        legend.key.size = unit(7,"points"),
        legend.position="top")+
  xlab("Assemblage")+
  guides(colour = guide_legend(nrow = 1, byrow = TRUE,
                               override.aes = list(alpha = 1)))


## Compare Functional richness with species richness

funric<-mymetrics %>% 
  mutate( zone=as.factor(case_when(clust==1~"A",
                                   clust==2~"B",
                                   clust==3~"C",
                                   clust==4~"D",
                                   clust==7~"E",
                                   clust==5~ "F",
                                   clust==6~ "G",
                                   clust==9~"H",
                                   clust==8~"I"))) %>% 
  mutate(zone=fct_relevel(zone,c("A","B","C","D","E","F","G","H","I")),
         Srich_st=Srich/max(Srich)) %>% 
  ggplot(aes(x=Srich,y=fric))+
  geom_point(aes(colour=zone))+
  scale_colour_manual(values=mypalette)+
  #geom_abline(slope=1,intercept=0)+
  theme_light()+
  theme(legend.title=element_blank(),
        legend.position=c(0.1,0.6),
        legend.background = element_blank())+
  xlab("Species richness")+ylab("Functional richness")
