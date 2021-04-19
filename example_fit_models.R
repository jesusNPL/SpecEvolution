
##### Usage #####
source("R/wrappers.R") # Set of functions that fit a set of evolutionary models using spectra

load("Data/oak_tree.rda") # Phylogeny
load("Data/oak_spec_jcb2016.rda") # Spectra data

oak_spec_jcb2016[1:5, 1:5]

View(oak_spec_jcb2016) # this data was obtained from Dudu's github

##### Data aggregation #####

# Given there are taxa with a different number of observations, let's aggregate the data
# This function return a list of four data.frames with the aggregated spectra by species
oak_spectra <- demon_AGG(spectra = oak_spec_jcb2016)

oak_spectra$Spec_mean[1:5, 1:5] # Mean spectra
oak_spectra$Spec_SE[1:5, 1:5] # Standard error spectra

##### Fit evolutionary models #####

# The aggregated data contains several columns, for practicality I'm selecting the first 10 bands
# You can also select the bands you want to evaluate by hand
# This function will return a list of three data.frames with all the parameters corresponding
# to the evolutionary model

# This is assuming NO measurement error
fit_spectra <- demon_Evol(spectra = oak_spectra$Spec_mean, # mean spectra data
                          tree = oak_tree, # phylogenetic tree
                          nBands = 10, # number of bands
                          NC = 2) # number of cores
fit_spectra

# This is assuming measurement error
fit_spectra_ME <- demon_Evol_ME(spectra = oak_spectra$Spec_mean, # mean spectra data
                                spec_ME = oak_spectra$Spec_SE, # SE spectra data
                                tree = oak_tree, # phylogenetic tree
                                nBands = 10, # number of bands
                                NC = 2) # number of cores
fit_spectra_ME

##### Compare the three evolutionary models using AIC ##### 

# This function returns a data.frame with the AIC weights for each evolutionary model
# that can be used to estimate the model evidence ratio

Comparison_spectra <- demon_ModSel(BM = fit_spectra_ME$BM_spectra$bm_AIC, 
                    OU = fit_spectra_ME$OU_spectra$ou_AIC, 
                    EB = fit_spectra_ME$EB_spectra$eb_AIC, 
                    nBands = 10, 
                    bandNames = fit_spectra_ME$BM_spectra$bandNames)
Comparison_spectra

##### Extract evidence ratio between models #####

# This function estimate the evidence ratio between the tree models based on the AIC weights
# The best model has NA value 

evidence_spectra <- demon_Evidence(Comparison = Comparison_spectra)
