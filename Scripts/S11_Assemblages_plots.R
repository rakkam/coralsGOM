##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 10: Make maps including information for coral assemblages
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Libraries

library(tidyverse)
library("ggsci")
library(raster)
library(sf)
library(ggpubr)
library(ggOceanMapsData)
library(ggOceanMaps)
library(ggnewscale)
library(oce)
library(formattable)

## Set directories
setwd("...")
localDir = "."
rawdataDir = file.path(localDir, "1_RawData")
modelDir = file.path(localDir, "2_Models")
resultDir = file.path(localDir, "results")
outputDir=file.path(localDir,"3_Outputs")
if (!dir.exists(resultDir)) dir.create(resultDir)

## Load data

load(paste(outputDir,"RCP_species_info.RData",sep="/"))

myclusters_final<-myclusters %>% 
  mutate( zone=as.factor(case_when(clust==1~"A",
                                   clust==2~"B",
                                   clust==3~"C",
                                   clust==4~"D",
                                   clust==7~"E",
                                   clust==5~ "F",
                                   clust==6~ "G",
                                   clust==9~"H",
                                   clust==8~"I"))) %>% 
  mutate(zone=fct_relevel(zone,c("A","B","C","D","E","F","G","H","I")))

## Make map

# Get crs and add depth

bath<-raster(file.path(rawdataDir,"misc/depth.grd"))

myclusters_sf<-myclusters_final%>% 
  st_as_sf(coords=c("Lon","Lat"),crs=crs(bath)) %>% 
  as_Spatial()

depth_clus<-extract(bath,myclusters_sf)

# Prepare basemap

load(file=paste(rawdataDir,"misc/shapefiles_for_basemap.RData",sep="/"))

mypalette<-c("#ABC4AB",
            "#855A5C",
            "#6E92A1",
            "#080E3A",
            "#E2C044",
            "#A87438",
            "#ABB557",
            "#27556C",
            "#8A8E91")
            

assem_map<-basemap(shapefiles=list(land = bs_land, glacier = NULL,bathy=bs_bathy),
        c(-73,-63, 38.5,45),bathymetry = TRUE,grid.col = NA,
        bathy.style = "poly_greys",land.col = "#eeeac4",bathy.alpha = 0.4,
        land.border.col = NA)+
  annotation_scale(location = "br",style="ticks",text_col="white",line_col="white") + 
  annotation_north_arrow(location = "tr", which_north = "true",height = unit(1, "cm"),
                         width = unit(1, "cm"),) +
  #  geom_point(data = model_expdes, aes(x = Lon, y = Lat), 
  #             color = "black",shape=21,size=2.5,fill="goldenrod2")+
  labs(x = NULL, y = NULL)+
  guides(fill = "none")+
  geom_spatial_point(data = myclusters_final, aes(x = Lon, y = Lat, color =zone),
                     crs=crs(bath),size=0.6)+
  #scale_color_uchicago(palette = c("light"))+
  #scale_color_npg(palette = "nrc")
  scale_colour_manual(values=mypalette)+
  theme(legend.position = "none",legend.title=element_blank())+
  guides(color=guide_legend(nrow=1, byrow=TRUE)) 

## Prepare T-S graph with assemblages
# read background temperature and salinity data for all the watercolumn

all_dat<-read.csv(paste(rawdataDir,"misc/all_wc_quadrat_data.csv",sep="/"))

muclust_coor<-myclusters_sf %>% 
  st_as_sf() %>% 
  st_transform(crs=4326)

st_coordinates(muclust_coor)

myclusters_clim<-cbind(myclusters_final,depth_clus) %>% 
  mutate(Lat=st_coordinates(muclust_coor)[,1],
         Lon=st_coordinates(muclust_coor)[,2],
           theta=swTheta(sal_mean,
    temperature = temp_mean,
    pressure = abs(depth_clus)),
    abs_sal=swAbsoluteSalinity(
      sal_mean,
      pressure = abs(depth_clus),
      latitude=Lon,longitude=Lat))


TS_graph<-all_dat %>% 
  ggplot(aes(x=sal,y=temp))+
  geom_point(colour="gray",size=0.2)+
  theme_light()+
  geom_point(data=myclusters_clim,aes(x=sal_mean,y=theta,colour=zone),size=0.1)+
  ylab("Temperature")+xlab("Salinity")+
  scale_colour_manual(values=mypalette)+
  scale_x_continuous(limits = c(32.7,35.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2,16), expand = c(0, 0))+
  theme(legend.title=element_text(size=7),
        legend.position="none",
        axis.title=element_text(size=7),
        axis.text=element_text(size=7))


TS_graph_blow<-all_dat %>% 
  ggplot(aes(x=sal,y=temp))+
  geom_point(colour="gray",size=0.2)+
  theme_light()+
  geom_point(data=myclusters_clim,aes(x=sal_mean,y=theta,colour=zone),size=0.1)+
  ylab(paste("Potential temperature (","\u00B0","C)",sep=""))+xlab("Absolute salinity")+
  scale_colour_manual(values=mypalette)+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=7),
        axis.text=element_text(size=7))+
  xlim(33.8,35.2)+ylim(2,9)

## Make table with average environmental variables for each assemblage

env_tab<-myclusters_clim%>%
  rename(Assemblage="zone") %>% 
  group_by(Assemblage) %>% 
  summarize(Depth=mean(depth_clus),
            Temperature=mean(temp_mean),
            Salinity=mean(sal_mean),
            Velocity=mean(vel_mean),
            Aspect=mean(aspect2000),
            Mud=mean(mud),
            TPI=mean(tpi20km)) %>%
  mutate(across(where(is.numeric),~round(.,digits=2))) %>% 
  #dplyr::arrange(desc(Depth)) %>% 
  rename(`Depth (m)`=Depth,
         `Temperature (C)`=Temperature,
         `Velocity (m/s)`=Velocity,
         `Mud (%)`=Mud)

mytab<-formattable(env_tab,align =rep("c",7),
            list(`Depth (m)`=color_tile("lightblue", "#d0867c"),
                 `Temperature (C)`=color_tile("lightblue", "#d0867c"),
                 `Salinity`=color_tile("lightblue", "#d0867c"),
                 `Velocity (m/s)`=color_tile("lightblue", "#d0867c"),
                 `Aspect`=color_tile("lightblue", "#d0867c"),
                 `Mud (%)`=color_tile("lightblue", "#d0867c"),
                 `TPI`=color_tile("lightblue", "#d0867c")))

## Prepare species heatmap

cluster_species_sum<-myclusters_final[,c(1:30,42)] %>% 
  group_by(zone) %>%
  summarize_all(mean) %>% 
  pivot_longer(cols=!zone,names_to="Species",values_to="prob_occ")

sp_heatmap<-ggplot(cluster_species_sum, aes(zone,Species)) +
  geom_tile(aes(fill = prob_occ),color = "gray87",linewidth=2)+
  geom_text(aes(label = round(prob_occ,2)), color = "black", size = 2)+
  scale_fill_gradient(low= "#fee6ce", high = "#e6550d",space = "Lab") +
  coord_fixed(ratio=0.4)+
  scale_x_discrete(position="top")+
  theme(axis.title=element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.title=element_blank(),
        legend.key.size = unit(1,"line"),
        legend.text=element_text(size=5))

