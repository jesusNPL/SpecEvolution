##### Usage #####
source("R/wrappers.R") # Set of functions that fit a set of evolutionary models using spectra

load("Data/oak_tree.rda") # Phylogeny
load("Data/oak_spec_jcb2016.rda") # Spectra data

# Given taht there are taxa with a different number of observations, let's aggregate the data
# This function return a list of four data.frames with the aggregated spectra by species
oak_spectra <- demon_AGG(spectra = oak_spec_jcb2016)

oak_spectra$Spec_mean[1:5, 1:5] # Mean spectra
oak_spectra$Spec_SE[1:5, 1:5] # Standard error spectra

dim(oak_spectra$Spec_mean)

## This is a function check the adequacy of traits (spectra) values in PCM.

ttt <- check_Adequacy(tree = oak_tree, 
                      spectra = oak_spectra$Spec_mean, 
                      nBands = 100, 
                      nSim = 1000)
