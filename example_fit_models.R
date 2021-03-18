
##### Usage #####
source("R/wrappers.R") # Set of functions that fit a set of evolutionary models using spectra

load("Data/oak_tree.rda") # Phylogeny
load("Data/oak_spec_jcb2016.rda") # Spectral data

View(oak_spec_jcb2016) # this data was obtained from Dudu's github

### Data aggregation
# Given there are taxa with a different number of observations, let's aggregate the data
# This function return a list of two data.frames with the aggregated spectra by species
oak_spectra <- demon_AGG(spectra = oak_spec_jcb2016)

### Fit evolutionary models 
# The aggregated data contains several columns, for simplicity I'm selecting the first 10 bands
# This function will return a list of three data.frames with all the parameters corresponding
# to the evolutionary model

fit_spectra <- demon_Evol(spectra = oak_spectra$Spec_mean, 
                          tree = oak_tree, 
                          nBands = 10)

### Compare the three evolutionary models using AIC 
# This function returns a data.frame with the AIC weights for each evolutionary model
# that can be used to estimate the model evidence ratio

Comparison_spectra <- demon_ModSel(BM = fit_spectra$BM_spectra$bm_AIC, 
                    OU = fit_spectra$OU_spectra$ou_AIC, 
                    EB = fit_spectra$EB_spectra$eb_AIC, 
                    nBands = 10, 
                    bandNames = fit_spectra$BM_spectra$bandNames)
