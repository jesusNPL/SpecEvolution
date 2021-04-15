
##### Auxiliary functions #####
SE <- function(x) {
  res <- sd(x)/sqrt(sum(!is.na(x)))
}

# Wrapper function that aggregate spectral data by taxa 
demon_AGG <- function(spectra) { 
  if ( ! ("plotrix" %in% installed.packages())) {install.packages("plotrix", dependencies = T)}
  
  Spectra_Mean <- aggregate( . ~ species, spectra, FUN = mean, na.rm = TRUE) 
  Spectra_SD <- aggregate( . ~ species, spectra, FUN = sd, na.rm = TRUE) 
  #Spectra_Median <- aggregate( . ~ species, FUN = spectra, median, na.rm = TRUE)
  Spectra_StdErr <- aggregate(. ~ species, spectra, FUN = SE) 
  # This is using the SE of the plotrix package
  Spectra_StdErr2 <- aggregate(. ~ species, spectra, FUN = plotrix::std.error)
  
  results <- list(Spec_mean = Spectra_Mean, Spec_SD = Spectra_SD, 
                  Spec_SE = Spectra_StdErr, Spec_SE2 = Spectra_StdErr2)
  return(results)  
}

# Wrapper to fit different evolutionary models

demon_Evol <- function(spectra, tree, nBands = 10, NC) {
  if ( ! ("geiger" %in% installed.packages())) {install.packages("geiger", dependencies = T)} 
  if ( ! ("svMisc" %in% installed.packages())) {install.packages("svMisc", dependencies = T)} 
  
  require(geiger)
  
  ## Inits
  nBands = nBands
  taxa <- spectra[, 1]
  spectrum <- spectra[, 2:(nBands + 1)]
  bandNames <- names(spectrum)[1:nBands]
  
  # White noise parameters
  
  # Brownian motion parameters
  bm_sigsq <- numeric(length = nBands)
  bm_z0 <- numeric(length = nBands)
  bm_lnL <- numeric(length = nBands)
  bm_nPar <- numeric(length = nBands)
  bm_AIC <- numeric(length = nBands)
  bm_AICc <- numeric(length = nBands)
  
  # Ornstein-Uhlenbeck parameters
  ou_alpha <- numeric(length = nBands)
  ou_sigsq <- numeric(length = nBands)
  ou_z0 <- numeric(length = nBands)
  ou_lnL <- numeric(length = nBands)
  ou_nPar <- numeric(length = nBands)
  ou_AIC <- numeric(length = nBands)
  ou_AICc <- numeric(length = nBands)
  ## Calculate the phylogenetic half-life -- 
  # how long does it take for half the information in the phylogeny to be erased
  halflife <- numeric(length = nBands)
  ## Compare this to total tree depth
  phy_halflife <- numeric(length = nBands)
  
  # Early burst parameters
  eb_a <- numeric(length = nBands)
  eb_sigsq <- numeric(length = nBands)
  eb_z0 <- numeric(length = nBands)
  eb_lnL <- numeric(length = nBands)
  eb_nPar <- numeric(length = nBands)
  eb_AIC <- numeric(length = nBands)
  eb_AICc <- numeric(length = nBands)
  
  for(i in 1:nBands) {
    svMisc::progress(i, max.value = nBands)
    
    band <- setNames(spectrum[, i], taxa)
      
    ## Brownian motion
    fit_bm <- geiger::fitContinuous(tree, band, model = "BM", ncores = NC)
    # Extract parameters
    bm_sigsq[i] <- fit_bm$opt$sigsq
    bm_z0[i] <- fit_bm$opt$z0
    bm_lnL[i] <- fit_bm$opt$lnL
    bm_nPar[i] <- fit_bm$opt$k
    bm_AIC[i] <- fit_bm$opt$aic
    bm_AICc[i] <- fit_bm$opt$aicc
    
    ## Ornstein-Uhlenbeck
    fit_ou <- geiger::fitContinuous(tree, band, model = "OU", ncores = NC)
    # Extract parameters
    ou_alpha[i] <- fit_ou$opt$alpha
    ou_sigsq[i] <- fit_ou$opt$sigsq
    ou_z0[i] <- fit_ou$opt$z0
    ou_lnL[i] <- fit_ou$opt$lnL
    ou_nPar[i] <- fit_ou$opt$k
    ou_AIC[i] <- fit_ou$opt$aic
    ou_AICc[i] <- fit_ou$opt$aicc
    
    ## Calculate the phylogenetic half-life -- 
    # how long does it take for half the information in the phylogeny to be erased
    halflife[i] <- log(2)/fit_ou$opt$alpha
    
    ## Compare this to total tree depth
    phy_halflife[i] <- halflife / max(branching.times(tree))
    
    ## Early burst
    fit_eb <- geiger::fitContinuous(tree, band, model = "EB", ncores = NC) 
    # Extract parameters
    eb_a[i] <- fit_eb$opt$a
    eb_sigsq[i] <- fit_eb$opt$sigsq
    eb_z0[i] <- fit_eb$opt$z0
    eb_lnL[i] <- fit_eb$opt$lnL
    eb_nPar[i] <- fit_eb$opt$k
    eb_AIC[i] <- fit_eb$opt$aic
    eb_AICc[i] <- fit_eb$opt$aicc
  }
  print(paste0("Three evolutionary models were fitted for ", nBands, 
               " bands assuming NO measurement error"))
  
  # Brownian motion parameters
  BM_spectra <- data.frame(bandNames, bm_sigsq, bm_z0, bm_lnL, 
                           bm_nPar, bm_AIC, bm_AICc)
  # Ornstein-Uhlenbeck parameters
  OU_spectra <- data.frame(bandNames, ou_alpha, ou_sigsq, ou_z0, ou_lnL, 
                           halflife, phy_halflife, 
                           ou_nPar, ou_AIC, ou_AICc)
  # Early burst parameters
  EB_spectra <- data.frame(bandNames, eb_a, eb_sigsq, eb_z0, eb_lnL, 
                           eb_nPar, eb_AIC, eb_AICc)
  
  results <- list(BM_spectra = BM_spectra, 
                  OU_spectra = OU_spectra, 
                  EB_spectra = EB_spectra)
  return(results)
}

# Wrapper to fit different evolutionary models considering measurement error

demon_Evol_ME <- function(spectra, spec_ME, tree, nBands = 10, NC) {
  if ( ! ("geiger" %in% installed.packages())) {install.packages("geiger", dependencies = T)} 
  if ( ! ("svMisc" %in% installed.packages())) {install.packages("svMisc", dependencies = T)} 
  
  require(geiger)
  
  ## Inits
  nBands = nBands
  taxa <- spectra[, 1]
  spectrum <- spectra[, 2:(nBands + 1)]
  bandNames <- names(spectrum)[1:nBands]
  
  # White noise parameters
  
  # Brownian motion parameters
  bm_sigsq <- numeric(length = nBands)
  bm_z0 <- numeric(length = nBands)
  bm_lnL <- numeric(length = nBands)
  bm_nPar <- numeric(length = nBands)
  bm_AIC <- numeric(length = nBands)
  bm_AICc <- numeric(length = nBands)
  
  # Ornstein-Uhlenbeck parameters
  ou_alpha <- numeric(length = nBands)
  ou_sigsq <- numeric(length = nBands)
  ou_z0 <- numeric(length = nBands)
  ou_lnL <- numeric(length = nBands)
  ou_nPar <- numeric(length = nBands)
  ou_AIC <- numeric(length = nBands)
  ou_AICc <- numeric(length = nBands)
  ## Calculate the phylogenetic half-life -- 
  # how long does it take for half the information in the phylogeny to be erased
  halflife <- numeric(length = nBands)
  ## Compare this to total tree depth
  phy_halflife <- numeric(length = nBands)
  
  # Early burst parameters
  eb_a <- numeric(length = nBands)
  eb_sigsq <- numeric(length = nBands)
  eb_z0 <- numeric(length = nBands)
  eb_lnL <- numeric(length = nBands)
  eb_nPar <- numeric(length = nBands)
  eb_AIC <- numeric(length = nBands)
  eb_AICc <- numeric(length = nBands)
  
  for(i in 1:nBands) {
    svMisc::progress(i, max.value = nBands)
    
    band <- setNames(spectrum[, i], taxa)
    ME <- setNames(spec_ME[, i], taxa)
    
    ## Brownian motion
    fit_bm <- geiger::fitContinuous(tree, band, model = "BM", SE = ME, ncores = NC)
    # Extract parameters
    bm_sigsq[i] <- fit_bm$opt$sigsq
    bm_z0[i] <- fit_bm$opt$z0
    bm_lnL[i] <- fit_bm$opt$lnL
    bm_nPar[i] <- fit_bm$opt$k
    bm_AIC[i] <- fit_bm$opt$aic
    bm_AICc[i] <- fit_bm$opt$aicc
    
    ## Ornstein-Uhlenbeck
    fit_ou <- geiger::fitContinuous(tree, band, model = "OU", SE = ME, ncores = NC)
    # Extract parameters
    ou_alpha[i] <- fit_ou$opt$alpha
    ou_sigsq[i] <- fit_ou$opt$sigsq
    ou_z0[i] <- fit_ou$opt$z0
    ou_lnL[i] <- fit_ou$opt$lnL
    ou_nPar[i] <- fit_ou$opt$k
    ou_AIC[i] <- fit_ou$opt$aic
    ou_AICc[i] <- fit_ou$opt$aicc
    
    ## Calculate the phylogenetic half-life -- 
    # how long does it take for half the information in the phylogeny to be erased
    halflife[i] <- log(2)/fit_ou$opt$alpha
    
    ## Compare this to total tree depth
    phy_halflife[i] <- halflife / max(branching.times(tree))
    
    ## Early burst
    fit_eb <- geiger::fitContinuous(tree, band, model = "EB", SE = ME, ncores = NC) 
    # Extract parameters
    eb_a[i] <- fit_eb$opt$a
    eb_sigsq[i] <- fit_eb$opt$sigsq
    eb_z0[i] <- fit_eb$opt$z0
    eb_lnL[i] <- fit_eb$opt$lnL
    eb_nPar[i] <- fit_eb$opt$k
    eb_AIC[i] <- fit_eb$opt$aic
    eb_AICc[i] <- fit_eb$opt$aicc
  }
  print(paste0("Three evolutionary models were fitted for ", nBands, 
               " bands assuming measurement error"))
  
  # Brownian motion parameters
  BM_spectra <- data.frame(bandNames, bm_sigsq, bm_z0, bm_lnL, 
                           bm_nPar, bm_AIC, bm_AICc)
  # Ornstein-Uhlenbeck parameters
  OU_spectra <- data.frame(bandNames, ou_alpha, ou_sigsq, ou_z0, ou_lnL, 
                           halflife, phy_halflife, 
                           ou_nPar, ou_AIC, ou_AICc)
  # Early burst parameters
  EB_spectra <- data.frame(bandNames, eb_a, eb_sigsq, eb_z0, eb_lnL, 
                           eb_nPar, eb_AIC, eb_AICc)
  
  results <- list(BM_spectra = BM_spectra, 
                  OU_spectra = OU_spectra, 
                  EB_spectra = EB_spectra)
  return(results)
}

# Wrapper to perform model selection
demon_ModSel <- function(BM, OU, EB, nBands = 10, bandNames) { 
  
  AIC_spectra <- list()
  
  for(k in 1:nBands) { 
    band <- bandNames[k]
    
    bm <- BM[k] 
    ou <- OU[k]
    eb <- EB[k]

    AICs <- c(bm, ou, eb)
    names(AICs) <- c("BM", "OU", "EB") 
    
    res <- geiger::aicw(AICs)
    
    res$Band <- band 
    res$Model <- c("BM", "OU", "EB")
    AIC_spectra[[k]] <- res
  }
  
  models <- do.call(rbind, AIC_spectra)
  
}

