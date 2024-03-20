# coralsGOM

This page includes the data and code for the manuscript:  
"Ocean circulation drives zonation of deep-water coral communities and their traits in the Northwest Atlantic"

The repository includes data and R scripts used to build [Hierarchical Models of Species Communities](https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc) to study the importance of environmental variables and coral traits on the assembly of coral communities in the NW Atlantic

## Data

The following datasets are included:

**4model_data_1km.Rdata**  
This is the curated dataset used to build the models. Includes the following vectors:
  - model_cordat: occurrences (presence/absence) of 30 coral genera in the study area
  - model_envdat: environmental variables in grid cells with coral occurrences
  - model_expdes: sampling design. Includes information for SurveyID, Dive, Lat, Lon for each grid cell with coral occurrences

**AxesTraits_1km.csv**  
Dataset created from the first Script S1. Includes the position of each coral genus in the three axes that describe the three-dimensional trait space

**Folder _misc_**: includes miscellaneous files that are needed for creating figures and maps
  - **Data_for_spatial_extrapolation.RData**  
Includes environmental variables, sampling effort and Lon/Lat for every grid cell in the   study area
  - **all_wc_quadrat_data.csv**  
Includes data of potential temperature and salinity for the whole water column above grid cells with coral occurrences. The data were obtained from the model fvcom 
  - **depth.grd** and **depth.gri**  
Raster files with depth data of all grid cells in the study area that are miscellaneously needed to construct spatial figures
  - **shapefiles_for_basemap.Rdata**  
Essential shapefiles to create maps. Acquired through the package XXX

## Scripts

The folder includes 12 script files. Each script produces data that are required for the scripts that follow, but the output of each script can be already found in the folder Data/3_Outputs. This allows to run each script independently if needed. The scripts to fit the models were modified from scripts available at the website of [Statistical Ecology group of the University of Helsinki](https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc)

**S1_Prepare_traits**  
Takes the dataset of coral traits as an input (included in *Data/1_RawData/4model_data_1km.RData*) and uses the package [mFD](https://cmlmagneville.github.io/mFD/) to create a three-dimensional trait space including all the coral traits. The script evaluates the quality of the trait space and estimates the optimum number of axes that can sufficiently describe it. Then it extracts the position of each of the 30 coral genera in the these trait axes (*Data/1_RawData/AxesTraits_1km.csv*), which will be used subsequently to construct the HMSC model.

**S2_define_models**  
Takes raw data (*Data/1_RawData/4model_data_1km.RData*) which include environmental variables, occurrences of 30 coral genera, coral traits and information on the experimental design and builds a HMSC model by using the package [Hmsc](https://github.com/hmsc-r/HMSC). It produces a file with the unfitted model (*Data/2_Models/unfitted_models.RData*)

**S3_fit_models**  
Fits the HMSC model by using Monte Carlo Chains  (4 chains, 250 samples, thinning by 100). Produces the fitted models (*Data/2_Models/models_thin_100_samples_250_chains_4.Rdata*). This code was run on a hypercomputer where it took 3 hours to complete, takes more than 48 hours on a conventional pc.

**S4_evaluate_convergence**  
Takes the fitted models (*Data/2_Models/models_thin_100_samples_250_chains_4.Rdata*) and produces graphs of the Gelman diagnostic to evaluate convergence of the different components of the models.

**S5_cross_validation**  
Takes the fitted models (*Data/2_Models/models_thin_100_samples_250_chains_4.Rdata*) and performs 10-fold cross validation. It uses parallel computing and it is preferably run on a supercomputer, as it takes very long time to run. Produces the file *Data/2_Models/MF_thin_100_samples_250_chains_4_nfolds_10.Rdata*

**S6_plot_model_results**  
Takes as input the fitted models (*Data/2_Models/models_thin_100_samples_250_chains_4.Rdata*), as well as the cross validation (*Data/2_Models/MF_thin_100_samples_250_chains_4_nfolds_10.Rdata*), and produces figures that summarize the model results. This includes figures of variance partitioning, explanatory and predicting power of the model, responses of each coral genus and coral traits to environmental variables.

**S7_plot_effects**  
Takes as input the fitted model (*Data/2_Models/models_thin_100_samples_250_chains_4.Rdata*) and produces figures for the total and marginal effects for each environmental variable, coral genera and traits.

**S8_spatial_extrapolation**  
Uses the package [dsmextra](https://github.com/densitymodelling/dsmextra) to evaluate the number of cells in the study area that have analogous conditions to those that were used to train the Hmsc model. Takes as input the environmental data for all the study area (*Data/1.RawData/misc/Data_for_spatial_extrapolation.RData*) and environmental data that were used to train the model (*Data/1.RawData/4model_data_1km.RData*). Produces a map with the locations of analogous and non-analogous areas, and a table with their percentages. It also produces a file with environmental data for all the analogous areas (*Data/3_Outputs/Analogue_areas_spatial_extrapolation.RData*), that is used later on to produce spatial predictions with the HMSC model.

**S9_produce_spatial_predictions**  
Uses the fitted HMSC model (*Data/2_Models/models_thin_100_samples_250_chains_4.Rdata*) to predict coral assemblages in all the analogous areas (*Data/3_Outputs/Analogue_areas_spatial_extrapolation.RData*). The script produces 1000 datasets and subsequently estimates average probability of coral occurrence, species richness and the standard deviation of species richness which is a measure of uncertainty, for each grid cell of the study area (analogous areas only). It also produces a map with species richness and its standard deviation.

**S10_Assemblages**  
The objective of this script is to group grid cells in coral assemblages based on their coral composition, and investigate the characteristics of each of these assemblages (such as species richness and trait diversity indices). The script takes a random sample of 5000 grid cells with analogous environments from the file (Analogue_areas_spatial_extrapolation.RData), and uses the HMSC model to produce 1000 datasets with predictions for coral occurrences for each grid cell of the study area. The average probability of coral occurrences is subsequently used as input to perform a Principal Component Analysis by using the packages [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html) and [FactoMiner](http://factominer.free.fr/). The script also uses the package [fundiversity](https://github.com/funecology/fundiversity) to calculate trait diversity indexes (functional richness, functional evenness, functional dispersion) for each of the 1000 datasets by using coral occurrence only (presence/absence), and subsequently estimates the average of each index for each grid cell. Produces a joint dataset with the assemblage number, the probability of occurrence of coral genera and trait index estimate for each grid cell (*Data/3_Outputs/RCP_species_info.RData*).

**S11_Assemblages_plots**  
Uses the dataset produced by the previous script (*Data/3_Outputs/RCP_species_info.RData*) to produce figures about the characteristics of each coral assemblage. It first creates a map with the study area and colours each grid cell to its respective assemblage. Subsequently, it creates a heatmap that shows the average occurrence of each coral genus in each assemblage. It also produces a salinity-temperature diagram that shows where each assemblage can be found in respect to the water column and the prevalent water masses in the area. Data for salinity and temperature were obtained by the model [fvcom-GOM](http://fvcom.smast.umassd.edu/research_projects/GB/gom.html). Lastly, it estimates the average of each environmental variable across the grid cells where each coral assemblage can be found.

**S12_Trait_graphs**  
This script uses the dataset (*Data/3_Outputs/RCP_species_info.RData*) to produce boxplots and dotplots showing the average Community Weighted Mean for each coral trait that was used in the study, as well as the species richness, functional richness, functional evenness and functional divergence of each assemblage.
