##----------------------------------------------------------------------------
## Deep-water coral communities in NW Atlantic
## Script 9: Performs Hierarchical Clustering to study assemblages
## (only using analogous areas)
## Author: Maria Rakka
##----------------------------------------------------------------------------

## Libraries

library(tidyverse)
library(foreach)
library(doParallel)
library(fundiversity)
library(Hmsc)
library(factoextra)
library(FactoMineR)
library(NbClust)

## Set and prepare directories
setwd("...")
localDir = "."
rawdataDir = file.path(localDir, "1_RawData")
modelDir = file.path(localDir, "2_Models")
outputDir=file.path(localDir,"3_Outputs")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

## Sample 5000 random grid cells from analogous areas
load(paste(outputDir,"Analogue_areas_spatial_extrapolation.RData",sep="/"))

analogue_dat<-analogue_dat %>% 
  rename(Lon=x,Lat=y)

set.seed(15)
an_slice<-slice_sample(analogue_dat,n=5000,weight_by=temp_mean)

## Load model and make predictions

load(paste(modelDir,"models_thin_100_samples_250_chains_4.RData"))

trial_env1<-an_slice %>% 
    dplyr::select(temp_mean,sal_mean,vel_mean,aspect2000,mud, tpi20km,sam_ef,ID=cell_nb) %>% 
    column_to_rownames("ID") %>% 
    as.data.frame()
  
xydat<-an_slice %>% 
    dplyr::select(Lon,Lat) %>% 
    as.matrix()
  
gradient_trial=prepareGradient(models$pres.abs.quad,
                                 XDataNew=trial_env1,
                                 sDataNew=list("ID"=xydat))
predict_trial=predict(models$pres.abs.quad,Gradient=gradient_trial,predictEtaMean = TRUE)
EpredY=Reduce("+",predict_trial)/length(predict_trial)

save(predict_trial,EpredY,file=file.path(outputDir,"Assemblage_prediction_outputs_raw.RData"))

## Perform PCA
# Attention the PCA takes long time to run
res.PCA<-PCA(EpredY,graph=FALSE)
ind <- get_pca_ind(res.PCA)

# Attention the NbClust function takes long time to run, may crash
res.nbclust <- NbClust(ind$coord, distance = "euclidean",
                       min.nc = 2, max.nc = 10, 
                       method = "complete", index ="all")

factoextra::fviz_nbclust(res.nbclust) + theme_minimal() +
  ggtitle("NbClust's optimal number of clusters")

save(res.nbclust,res.PCA,file=paste(resultDir,"PCA_results.RData",sep="/"))

res.nbclust$Best.nc

mypalette_den<-c("#ABC4AB",
                          "#855A5C",
                          "#6E92A1",
                          "#080E3A",
                          "#E2C044",
                          "#A87438",
                          "#ABB557",
                          "#27556C",
                          "#8A8E91")
                          
res.HCPC<-HCPC(res.PCA,nb.clust=9,consol=FALSE,graph=FALSE)

fviz_pca_var(res.PCA, col.var="contrib",select.var = list(contrib=10),
             gradient.cols = c("#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

fviz_cluster(res.HCPC,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             geom = "point",
             palette = "npg",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)


plot.HCPC(res.HCPC,choice='tree',title='Hierarchical tree')
plot(res.HCPC, choice='tree',cex = 0.6)


# Manual dendrogram, ATTENTION takes long time to run

fviz_dend2 <-function (x, k = NULL, h = NULL, k_colors = NULL, palette = NULL, 
            show_labels = TRUE, color_labels_by_k = TRUE, label_cols = NULL, 
            labels_track_height = NULL, repel = FALSE, lwd = 0.7, type = c("rectangle", 
                                                                           "circular", "phylogenic"), phylo_layout = "layout.auto", 
            rect = FALSE, rect_border = "gray", rect_lty = 2, rect_fill = FALSE, 
            lower_rect, horiz = FALSE, cex = 0.8, main = "Cluster Dendrogram", 
            xlab = "", ylab = "Height", sub = NULL, ggtheme = theme_classic(), 
            ...) {
    if (missing(k_colors) & !is.null(palette)) {
      k_colors <- palette
      palette <- NULL
    }
    if (!color_labels_by_k & is.null(label_cols)) 
      label_cols <- "black"
        type <- match.arg(type)
        circular <- type == "circular"
        phylogenic <- type == "phylogenic"
        rectangle <- type == "rectangle"
        if (inherits(x, "HCPC")) {
          k <- length(unique(x$data.clust$clust))
          x <- x$call$t$tree
        }
        if (inherits(x, "hcut")) {
          k <- x$nbclust
          dend <- as.dendrogram(x)
          method <- x$method
        }
        else if (inherits(x, "hkmeans")) {
          k <- length(unique(x$cluster))
          dend <- as.dendrogram(x$hclust)
          method <- x$hclust$method
        }
        else if (inherits(x, c("hclust", "agnes", "diana"))) {
          dend <- as.dendrogram(x)
          method <- x$method
        }
        else if (inherits(x, "dendrogram")) {
          dend <- x
          method <- ""
        }
        else stop("Can't handle an object of class ", paste(class(x), 
                                                            collapse = ", "))
        if (is.null(method)) 
          method <- ""
        else if (is.na(method)) 
          method <- ""
        if (is.null(sub) & method != "") 
          sub = paste0("Method: ", method)
        if (!is.null(dendextend::labels_cex(dend))) 
          cex <- dendextend::labels_cex(dend)
        dend <- dendextend::set(dend, "labels_cex", cex)
        dend <- dendextend::set(dend, "branches_lwd", lwd)
        k <- factoextra:::.get_k(dend, k, h)
        if (!is.null(k)) {
          if (ggpubr:::.is_col_palette(k_colors)) 
            k_colors <- ggpubr:::.get_pal(k_colors, k = k)
          else if (is.null(k_colors)) 
            k_colors <- ggpubr:::.get_pal("default", k = k)
          dend <- dendextend::set(dend, what = "branches_k_color", 
                                  k = k, value = k_colors)
          if (color_labels_by_k) 
            dend <- dendextend::set(dend, "labels_col", k = k, 
                                    value = k_colors)
        }
        if (!is.null(label_cols)) {
          dend <- dendextend::set(dend, "labels_col", label_cols)
        }
        leaflab <- ifelse(show_labels, "perpendicular", "none")
        if (xlab == "") 
          xlab <- NULL
        if (ylab == "") 
          ylab <- NULL
        max_height <- max(dendextend::get_branches_heights(dend))
        if (missing(labels_track_height)) 
          labels_track_height <- max_height/8
        if (max_height < 1) 
          offset_labels <- -max_height/100
        else offset_labels <- -0.1
        if (rectangle | circular) {
          p <- factoextra:::.ggplot_dend(dend, type = "rectangle", offset_labels = offset_labels, 
                                         nodes = FALSE, ggtheme = ggtheme, horiz = horiz, 
                                         circular = circular, palette = palette, labels = show_labels, 
                                         label_cols = label_cols, labels_track_height = labels_track_height, 
                                         ...)
          if (!circular) 
            p <- p + labs(title = main, x = xlab, y = ylab)
        }
        else if (phylogenic) {
          p <- .phylogenic_tree(dend, labels = show_labels, label_cols = label_cols, 
                                palette = palette, repel = repel, ggtheme = ggtheme, 
                                phylo_layout = phylo_layout, ...)
        }
        if (circular | phylogenic | is.null(k)) 
          rect <- FALSE
        if (rect_fill & missing(rect_lty)) 
          rect_lty = "blank"
        if (missing(lower_rect)) 
          lower_rect = -(labels_track_height + 0.5)
        if (rect) {
          p <- p + rect_dendrogram(dend, k = k, palette = rect_border, 
                                   rect_fill = rect_fill, rect_lty = rect_lty, size = lwd, 
                                   lower_rect = lower_rect)
        }
        attr(p, "dendrogram") <- dend
        structure(p, class = c(class(p), "fviz_dend"))
        return(p)
  }

rect_dendrogram <- function(dend, k = NULL, h = NULL, k_colors = NULL, palette = NULL, 
                            rect_fill = FALSE, rect_lty = 2, lower_rect = -1.5, ...) {
  if (missing(k_colors) & !is.null(palette)) 
    k_colors <- palette
  prop_k_height <- 0.5
  if (!dendextend::is.dendrogram(dend)) 
    stop("x is not a dendrogram object.")
  k <- factoextra:::.get_k(dend, k, h)
  tree_heights <- dendextend::heights_per_k.dendrogram(dend)[-1]
  tree_order <- stats::order.dendrogram(dend)
  if (is.null(k)) 
    stop("specify k")
  if (k < 2) {
    stop(gettextf("k must be between 2 and %d", length(tree_heights)), 
         domain = NA)
  }
  cluster <- dendextend::cutree(dend, k = k)
  clustab <- table(cluster)[unique(cluster[tree_order])]
  m <- c(0, cumsum(clustab))
  which <- 1L:k
  xleft <- ybottom <- xright <- ytop <- list()
  for (n in seq_along(which)) {
    next_k_height <- tree_heights[names(tree_heights) == 
                                    k + 1]
    if (length(next_k_height) == 0) {
      next_k_height <- 0
      prop_k_height <- 1
    }
    xleft[[n]] = m[which[n]] + 0.66
    ybottom[[n]] = lower_rect
    xright[[n]] = m[which[n] + 1] + 0.33
    ytop[[n]] <- tree_heights[names(tree_heights) == k] * 
      prop_k_height + next_k_height * (1 - prop_k_height)
  }
  df <- data.frame(xmin = unlist(xleft), ymin = unlist(ybottom), 
                   xmax = unlist(xright), ymax = unlist(ytop), stringsAsFactors = TRUE)
  color <- k_colors
  if (all(color == "cluster"))
    color <- "default"
  if (ggpubr:::.is_col_palette(color)) 
    color <- ggpubr:::.get_pal(color, k = k)
  else if (length(color) > 1 & length(color) < k) {
    color <- rep(color, k)[1:k]
  }
  if (rect_fill) {
    fill <- color
    alpha <- 0.2
  }
  else {
    fill <- "transparent"
    alpha <- 0
  }
  df$color <- color
  df$cluster <- as.factor(paste0("c", 1:k))
  ggpubr::geom_exec(geom_rect, data = df, xmin = "xmin", ymin = "ymin", 
                    xmax = "xmax", ymax = "ymax", fill = fill, color = color, 
                    linetype = rect_lty, alpha = alpha, ...)
}


dendrogram_graph<-fviz_dend2(res.HCPC, k = 9, 
           cex = 0.6, # label size
           k_colors = mypalette_den,
           # color_labels_by_k = TRUE, # color labels by groups
           rect = TRUE, # Add rectangle around groups
           rect_border = mypalette_den,
           rect_fill = TRUE)

ggsave(dendrogram_graph,file=paste(resultDir,"dendrogram.png",sep="/"))

## PCA biplot

mypalette<-c("#ABC4AB",
                      "#855A5C",
                      "#6E92A1",
                      "#080E3A",
                      "#E2C044",
                      "#A87438",
                      "#ABB557",
                      "#27556C",
                      "#8A8E91")


allcoords<-as.data.frame(ind$coord)
allcoords$clust<-res.HCPC$data.clust$clust

allcoords<-allcoords %>% 
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


ggplot(allcoords,aes(x=Dim.1,y=Dim.2,colour=zone))+
  geom_point()+
  scale_colour_manual(values=mypalette,name="Assemblage")+
  xlim(-8,8)+
  theme_light()+
  xlab("Dim 1 (48%)")+ylab("Dim 2 (22.5%)")+
  theme(legend.position=c(0.1,0.6),
        legend.background=element_rect(colour="gray"))


#--------------------------------------------------------------------------
## Calculate trait indices
#--------------------------------------------------------------------------

mytraits<-models$pres.abs.quad$Tr[,c(2:4)]


mytrait_fun<-function(x){
  fric<-fd_fric(mytraits, x,stand = TRUE) 
  feve<-fd_feve(mytraits, x)
  fdis<-fd_fdis(mytraits, x)
  fdiv<-fd_fdiv(mytraits, x)
  return(cbind(fric=fric$FRic,feve=feve$FEve,fdis=fdis$FDis,fdiv=fdiv))
}


## Run each metric in parallel separately

mytrait_fun(predict_trial[[1]])

trial1<-predict_trial[1:250]
trial2<-predict_trial[251:500]
trial3<-predict_trial[501:750]
trial4<-predict_trial[751:1000]

cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

rf <- foreach(x=iter(trial1,chunksize=5), .packages=c("fundiversity")) %dopar%
  mytrait_fun(x)
stopCluster(cl)

cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

rf2 <- foreach(x=iter(trial2,chunksize=5), .packages=c("fundiversity")) %dopar%
  mytrait_fun(x)
stopCluster(cl)

cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
rf3 <- foreach(x=iter(trial3,chunksize=5), .packages=c("fundiversity")) %dopar%
  mytrait_fun(x)
stopCluster(cl)

cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
rf4 <- foreach(x=iter(trial4,chunksize=5), .packages=c("fundiversity")) %dopar%
  mytrait_fun(x)
stopCluster(cl)

all_rfs<-c(rf,rf2,rf3,rf4) 

all_rfs_cor<-rapply( all_rfs, f=function(x) ifelse(is.na(x),0,x), how="replace" ) %>% 
  map(~dplyr::select(.,!fdiv.site)) %>%
  map(~rename(.,fdiv=fdiv.FDiv)) 

trait_indices<-  Reduce("+",all_rfs_cor)/length(all_rfs_cor)


#--------------------------------------------------------------------------
## Other metrics (this was done before for all area, but now redoing for the 
## 5000 sampled grid cells)
#--------------------------------------------------------------------------

getS = function(p){return(rowSums(p))}
aS = simplify2array(lapply(X = predict_trial, FUN = getS))
dim(aS)
Srich = apply(aS, 1, mean)
sdS = sqrt(apply(aS, 1, var))
S5 = apply(aS > 5, 1, mean)
S10 = apply(aS > 10, 1, mean)
CWM = (EpredY %*% models$pres.abs.quad$Tr)/matrix(rep(Srich, models$pres.abs.quad$nt), ncol = models$pres.abs.quad$nt)

mymetrics<-cbind(data.frame(Srich,sdS,S5,S10,CWM),clust=myclusters$clust,trait_indices)

save(myclusters,mymetrics,file=file.path(outputDir,"RCP_species_info.RData"))

