##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 6: Plot model results and variance partitioning
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Libraries
library(Hmsc)
library(tidyverse)
library(ggpubr)
library(ggsci)

## Set directory

setwd("...")
localDir = "."
modelDir = file.path(localDir, "2_Models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

## Load fitted models
load(paste(modelDir,"models_thin_100_samples_250_chains_4.RData",sep="/"))

## Plot model results (environmental variables and traits)

postBeta<-getPostEstimate(models$pres.abs, parName="Beta")
postBeta [2:3]%>% 
  map(.,~as.data.frame(.)) %>% 
  map(.,~mutate(.,across(everything(),~na_if(.,.>0.95))))


givemeneg<-function(x){
  newx<-ifelse(is.na(x),NA,x*(-1))
  return(newx)
}

varnames<-c("Intercept","Salinity","Temperature",
            "Temperature^2","Velocity","TPI","Mud","Aspect","SE")
postbeta_df<-postBeta [2:3]%>% 
  map(.,~as.data.frame(.)) %>% 
  map(.,~mutate(.,across(everything(), ~ifelse( .<0.95, NA, .)))) %>% 
  list_rbind() %>% 
  mutate(supp=rep(c("pos","neg"),each=9))%>% 
  #mutate(Stauropathes=ifelse(supp=="neg",givemeneg(Stauropathes),Stauropathes))
  mutate(across(where(is.numeric),~ifelse(supp=="neg",givemeneg(.),.)),
         myvars=fct_relevel(as_factor(rep(varnames,2)),varnames)) %>%
  #select(!supp) %>% 
  pivot_longer(cols=!c("myvars","supp"),names_to="Species",values_to="Support") %>% 
  group_by(Species,myvars) %>% 
  summarize(all_sup=sum(Support,na.rm=TRUE)) %>%
  drop_na() %>% 
  mutate(all_sup=na_if(all_sup,0))


betas_graph<-ggplot(postbeta_df, aes(x=myvars,y=Species)) +
  geom_tile(aes(fill = all_sup))+
  #geom_text(aes(label =round(all_sup,2)), color = "gray", size = 3)+
  scale_fill_gradient2(low = '#006699',
                       high = '#b93a26',na.value = "white",breaks=c(-1,0,1)) +
  #scale_fill_gradient(low = "blue", high = "red",limits = c(-1, 1))+
  coord_fixed(ratio=0.3)+
  #scale_x_discrete(position="top")+
  theme(axis.title=element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=6),
        axis.text.y=element_text(colour="black",size=6),
        legend.title=element_blank(),
        legend.key.size = unit(9,"points"))+
  scale_y_discrete(limits=rev)+
  scale_x_discrete(labels=c("Intercept","Salinity","Temperature",
                             bquote(Temperature^2),"Velocity",
                            "TPI","Mud","Aspect","SE"))

postGamma<-getPostEstimate(models$pres.abs, parName = "Gamma")

postGamma_df<-postGamma [2:3]%>% 
  map(.,~as.data.frame(.)) %>% 
  map(.,~mutate(.,across(everything(), ~ifelse( .<0.95, NA, .)))) %>% 
  list_rbind() %>% 
  mutate(supp=rep(c("pos","neg"),each=9))%>% 
  #mutate(Stauropathes=ifelse(supp=="neg",givemeneg(Stauropathes),Stauropathes))
  mutate(across(where(is.numeric),~ifelse(supp=="neg",givemeneg(.),.)),
         myvars=fct_relevel(as_factor(rep(varnames,2)),varnames)) %>%
  #select(!supp) %>% 
  pivot_longer(cols=!c("myvars","supp"),names_to="Species",values_to="Support") %>% 
  group_by(Species,myvars) %>% 
  summarize(all_sup=sum(Support,na.rm=TRUE)) %>%
  drop_na() %>% 
  mutate(all_sup=na_if(all_sup,0))

traits_graph<-ggplot(postGamma_df, aes(x=myvars,y=Species)) +
  geom_tile(aes(fill = all_sup))+
  #geom_text(aes(label =round(all_sup,2)), color = "gray", size = 3)+
  scale_fill_gradient2(low = '#006699',
                       high = '#b93a26',na.value = "white",breaks=c(-1,0,1)) +
  #scale_fill_gradient(low = "blue", high = "red",limits = c(-1, 1))+
  coord_fixed(ratio=0.7)+
  #scale_x_discrete(position="top")+
  theme(axis.title=element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=6),
        axis.text.y=element_text(colour="black",size=6),
        legend.title=element_blank(),
        legend.key.size = unit(9,"points"))+
  scale_y_discrete(labels=c("Intercept","PC1","PC2","PC3"))+
  scale_x_discrete(labels=c("Intercept","Salinity","Temperature",
                            bquote(Temperature^2),"Velocity",
                            "TPI","Mud","Aspect","SE"))

ggarrange(betas_graph,traits_graph,labels=c("A","B"), ncol=1,
          heights = c(4.2, 2),align="v",common.legend=TRUE,legend="bottom")

##Variance partitioning

VP = computeVariancePartitioning(models$pres.abs)
VP$R2T
VP$R2T$Y
newvardat<-
as.data.frame(VP$vals) %>% 
  rownames_to_column("Variable") %>% 
  mutate(Variable=fct_recode(Variable,"Salinity"="sal_mean",
                             "Temperature"="poly(temp_mean, degree = 2, raw = TRUE)",
                             "Velocity"="vel_mean",
                             "TPI"="tpi20km",
                             "Mud"="mud",
                             "Aspect"="aspect2000",
                             "Sam.effort"="sam_ef")) %>% 
  pivot_longer(cols=!Variable,names_to = "Genus",values_to = "varpar") %>% 
  mutate(Variable=fct_relevel(Variable,"Salinity","Temperature","Velocity",
                              "TPI","Aspect","Mud"),
         varpar=100*varpar)

round(100 * rowMeans(VP$vals), 1)

newvardat %>% 
  ggplot(aes(x=Genus,y=varpar))+
  geom_bar(aes(fill=Variable),stat="identity",position="stack")+
  scale_fill_simpsons(labels=c("Salinity (37.0)",
                              "Temperature (25.3)",
                              "Velocity (4.5)",
                              "TPI (4.7)",
                              "Aspect (3.5)",
                              "Mud (3.3)",
                              "Random: Dive (6.2)",
                              "Random ID (8.8)",
                              "Sam. effort (6.6)"))+
  theme_minimal()+
  ylab("Variance (%)")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(),
        legend.title=element_blank())

head(models$pres.abs$X)
VP2 = computeVariancePartitioning(models$pres.abs,group=c(1,1,1,1,1,2,2,2,3),
                                  groupnames = c("Oceanographic","Terrain","Sam.effort"))

newvardat2<-
  as.data.frame(VP2$vals) %>% 
  rownames_to_column("Variable") %>% 
  pivot_longer(cols=!Variable,names_to = "Genus",values_to = "varpar") %>% 
  mutate(varpar=100*varpar,
         Variable=fct_relevel(Variable,"Oceanographic","Terrain"))

round(100 * rowMeans(VP2$vals), 1)

newvardat2 %>% 
  ggplot(aes(x=Genus,y=varpar))+
  geom_bar(aes(fill=Variable),stat="identity",position="stack")+
  scale_fill_simpsons(labels=c("Oceanographic (69.8)",
                               "Terrain (9.6)",
                               "Random: Dive (6.2)",
                               "Random ID (10.3)",
                               "Sam. effort (5.5)"))+
  theme_minimal()+
  ylab("Variance (%)")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(),
        legend.title=element_blank())

##Explanatory power (R2)

preds = computePredictedValues(models$pres.abs)
MF = evaluateModelFit(hM=models$pres.abs, predY=preds)

r2_gen<-data.frame(R2=MF$TjurR2,Genus=models$pres.abs$spNames)

ggplot(r2_gen,aes(y=reorder(Genus,R2),x=R2))+
  geom_bar(stat="identity")+
  theme_light()+
  theme(axis.title.y=element_blank())+
  xlab(bquote('Explanatory '~R^2))

##predicting power

#load data from cross-validation
load(paste(modelDir,"MF_thin_100_samples_250_chains_4_nfolds_10.RData",sep="/"))

#MFCV predicting power

all_modelfit<-data.frame(Species=models$pres.abs$spNames,
                         exR=MF$pres.abs$RMSE,exAUC=MF$pres.abs$AUC,exTjurR2=MF$pres.abs$TjurR2,
                         predR=MFCV$pres.abs$RMSE,predAUC=MFCV$pres.abs$AUC,predTjurR2=MFCV$pres.abs$TjurR2,
                         SpOcc=colSums(data.frame(models$pres.abs$Y) ))

pred_tjur<-ggplot(all_modelfit,aes(x=predTjurR2,y=reorder(Species,predTjurR2)))+
  geom_bar(stat="identity")+
  #geom_vline(xintercept=0.5, linetype="dashed", color = "red")+
  theme_light()+
  ylab("")+xlab(bquote('Predictive'~R^2))

pred_auc<-ggplot(all_modelfit,aes(x=predAUC,y=reorder(Species,predAUC)))+
  geom_bar(stat="identity")+
  #geom_vline(xintercept=0.5, linetype="dashed", color = "red")+
  theme_light()+
  ylab("")+xlab("Predictive AUC")

all_modelfit %>% 
  arrange(SpOcc)
       