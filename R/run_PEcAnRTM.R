##### Function that prepare data for being used in PEcAnRTM #####
# spec = spectra data

prepareSpec_PEcAnRTM <- function(spec, saveSpec = FALSE, 
                                 path = "Data/", 
                                 name = "OH") {
  library(spectrolab)
  # Guess splice
  spec_splice <- guess_splice_at(spec)
  # Match sensors
  spec_matched <- match_sensors(x = spec, 
                                splice_at = spec_splice, 
                                interpolate_wvl = c(5, 1)) 
  # Resample spectra to 1 nm
  spec <- spectrolab::resample(spec_matched, new_bands = seq(400, 2500, 1)) 
  spec_dta <- as.data.frame(spec, fix_names = "none", metadata = FALSE)
  
  if(saveSpec == TRUE)
  {
    save(spec, file = paste0(path, "spec_", name, ".RData"))
  }
  
  dta <- list(spec, spec_dta)
  
  return(dta)
}

# Usage
#specOH <- prepareSpec_PEcAnRTM(spec = bbb[[1]], saveSpec = FALSE)

#dt <- specOH[[2]]

##### Function to run PEcAnRTM #####

# specData = spectra data in data.frame format
# specRange = if range of spectra values, i.e., columns that have spectra values
# sppNames = vector of species names

runPEcAnRTM <- function(specData, specRange = c(2:2102), sppNames, 
                        nChains = 4, iters = 10000) { 
  library(PEcAnRTM)
  ### Inits
  params <- c(
    "N" = 1.4, "Cab" = 40, "Car" = 15,
    "Cbrown" = 0.5, "Cw" = 0.002, "Cm" = 0.004
  )
  
  ### PEcAnRTM paramenters
  invert.options <- default.settings.prospect
  invert.options$model
  
  invert.options$n.tries <- 1      # Number of attempts
  invert.options$nchains <- nChains      # Number of MCMC chains
  invert.options$ngibbs <- iters    # Number of iterations per chain
  invert.options$burnin <- iters/5    # Length of burnin period
  invert.options$do.lsq.first <- TRUE # Initialize with results from a fast least-squares optimization algorithm
  
  # Extra params
  wls <- 400:2500
  
  waves <- paste0("Wave_", wls)
  
  traitMeans <- list()
  traitMedians <- list()
  
  for(i in 1:length(sppNames)) {
    
    specTMP <- as.numeric(specData[i, ][specRange])
    names(specTMP) <- waves
    
    inversion.output <- PEcAnRTM::invert.auto(observed = specTMP,
                                              invert.options = invert.options,
                                              quiet = TRUE) 
    # Extract mean values
    means.tmp <- unlist(inversion.output$results[grep("mu", 
                                                     names(inversion.output$results))])[1:5]
    
    means.tmp <- as.data.frame(means.tmp)
    means.tmp$Species <- sppNames[i]
    
    traitMeans[[i]] <- means.tmp
    # Extract median values
    median.tmp <- unlist(inversion.tmp$results[grep("med", 
                                                      names(inversion.tmp$results))])[1:5]
    
    median.tmp <- as.data.frame(median.tmp)
    median.tmp$Species <- sppNames[i]
    
    traitMedians[[i]] <- median.tmp
  }
  ## Save trait means
  traitMeans <- do.call(rbind, traitMeans)
  ## Save trait medians
  traitMedians <- do.call(rbind, traitMedians)
  
  results <- list(traitMeans, traitMedians)
  
  return(results)
  
}

##### Usage

#load("Data/spec_OH.RData")
dta4 <- as.data.frame(spec, fix_names = "none", metadata = FALSE)

oaks <- dta4$sample_name

traitOH <- runPEcAnRTM(specData = dta4, specRange = c(2:2102), 
                       sppNames = oaks, nChains = 2, iters = 5000)



