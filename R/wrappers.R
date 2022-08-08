
##### Auxiliary functions #####
SE <- function(x) {
  res <- sd(x) / sqrt(sum(!is.na(x)))
}

# Wrapper function that aggregate spectral data by taxa
demon_AGG <- function(spectra) {
  if (!("plotrix" %in% installed.packages())) {
    install.packages("plotrix", dependencies = T)
  }

  Spectra_Mean <- aggregate(. ~ species, spectra, FUN = mean, na.rm = TRUE)
  Spectra_SD <- aggregate(. ~ species, spectra, FUN = sd, na.rm = TRUE)
  # Spectra_Median <- aggregate( . ~ species, FUN = spectra, median, na.rm = TRUE)
  Spectra_StdErr <- aggregate(. ~ species, spectra, FUN = SE)
  # This is using the SE of the plotrix package
  Spectra_StdErr2 <- aggregate(. ~ species, spectra, FUN = plotrix::std.error)

  results <- list(
    Spec_mean = Spectra_Mean, Spec_SD = Spectra_SD,
    Spec_SE = Spectra_StdErr, Spec_SE2 = Spectra_StdErr2
  )
  return(results)
}

##### Check adequacy #####

## This is a function that check the adequacy of traits (spectra) values in PCM.
# If the data is adequate it don't differ from the null distribution.

check_Adequacy <- function(tree, spectra, nBands = 100, nSim = 1000) {
  require(arbutus)
  require(geiger)

  ## Inits
  nBands <- nBands
  taxa <- spectra[, 1]
  spectrum <- spectra[, 2:(nBands + 1)]
  bandNames <- names(spectrum)[1:nBands]

  obs_tbl <- list()
  p_vals <- list()

  for (i in 1:nBands) {
    svMisc::progress(i, max.value = nBands)

    band <- setNames(spectrum[, i], taxa)

    treeData <- treedata(tree, band, sort = TRUE, warnings = FALSE)
    # This is annoying, tree and data should be in the same order
    band_order <- band[match(tree$tip.label, names(band))]

    unity <- make_unit_tree(tree, data = band_order)

    ## calculate default test stats on observed data
    obs <- calculate_pic_stat(unity, stats = NULL)
    # obs_tmp <- obs
    # obs_tmp$Band <- bandNames[i]
    obs_tbl[[i]] <- obs
    ## simulate data on unit.tree
    sim.dat <- simulate_char_unit(unity, nsim = nSim)
    ## calculate default test stats on simulated data
    sims <- calculate_pic_stat(sim.dat, stats = NULL)
    ## compare simulated to observed test statistics
    res <- compare_pic_stat(obs, sims)
    p_vals[[i]] <- res$p.values
  }

  obs_tbl <- do.call(rbind, obs_tbl)
  obs_tbl$Band <- bandNames
  p_vals <- data.frame(do.call(rbind, p_vals))
  p_vals$Band <- bandNames

  results <- list(adequacy = obs_tbl, pVals = p_vals)

  return(results)
}

# Wrapper to fit different evolutionary models

demon_Evol <- function(spectra, tree, nBands = 10, NC) {
  if (!("geiger" %in% installed.packages())) {
    install.packages("geiger", dependencies = T)
  }
  if (!("svMisc" %in% installed.packages())) {
    install.packages("svMisc", dependencies = T)
  }

  require(geiger)

  ## Inits
  nBands <- nBands
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

  for (i in 1:nBands) {
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
    halflife[i] <- log(2) / fit_ou$opt$alpha

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
  print(paste0(
    "Three evolutionary models were fitted for ", nBands,
    " bands assuming NO measurement error"
  ))

  # Brownian motion parameters
  BM_spectra <- data.frame(
    bandNames, bm_sigsq, bm_z0, bm_lnL,
    bm_nPar, bm_AIC, bm_AICc
  )
  # Ornstein-Uhlenbeck parameters
  OU_spectra <- data.frame(
    bandNames, ou_alpha, ou_sigsq, ou_z0, ou_lnL,
    halflife, phy_halflife,
    ou_nPar, ou_AIC, ou_AICc
  )
  # Early burst parameters
  EB_spectra <- data.frame(
    bandNames, eb_a, eb_sigsq, eb_z0, eb_lnL,
    eb_nPar, eb_AIC, eb_AICc
  )

  results <- list(
    BM_spectra = BM_spectra,
    OU_spectra = OU_spectra,
    EB_spectra = EB_spectra
  )
  return(results)
}

# Wrapper to fit different evolutionary models considering measurement error

demon_Evol_ME <- function(spectra, spec_ME, tree, nBands = 10, NC) {
  if (!("geiger" %in% installed.packages())) {
    install.packages("geiger", dependencies = T)
  }
  if (!("svMisc" %in% installed.packages())) {
    install.packages("svMisc", dependencies = T)
  }

  require(geiger)

  ## Inits
  nBands <- nBands
  taxa <- spectra[, 1]
  spectrum <- spectra[, 2:(nBands + 1)]
  spec_ME <- spec_ME[, 2:(nBands + 1)]

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

  for (i in 1:nBands) {
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
    halflife[i] <- log(2) / fit_ou$opt$alpha

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
  print(paste0(
    "Three evolutionary models were fitted for ", nBands,
    " bands assuming measurement error"
  ))

  # Brownian motion parameters
  BM_spectra <- data.frame(
    bandNames, bm_sigsq, bm_z0, bm_lnL,
    bm_nPar, bm_AIC, bm_AICc
  )
  # Ornstein-Uhlenbeck parameters
  OU_spectra <- data.frame(
    bandNames, ou_alpha, ou_sigsq, ou_z0, ou_lnL,
    halflife, phy_halflife,
    ou_nPar, ou_AIC, ou_AICc
  )
  # Early burst parameters
  EB_spectra <- data.frame(
    bandNames, eb_a, eb_sigsq, eb_z0, eb_lnL,
    eb_nPar, eb_AIC, eb_AICc
  )

  results <- list(
    BM_spectra = BM_spectra,
    OU_spectra = OU_spectra,
    EB_spectra = EB_spectra
  )
  return(results)
}

# Wrapper to perform model selection
demon_ModSel <- function(BM, OU, EB, nBands = 10, bandNames) {
  AIC_spectra <- list()

  for (k in 1:nBands) {
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

  models <- data.frame(do.call(rbind, AIC_spectra))
  names(models) <- c("AIC", "deltaAIC", "AICw", "Band", "Model")
  rownames(models) <- NULL

  models <- models[, c(5, 4, 1, 2, 3)]

  return(models)
}

# Function to calculate the Evidence ratio

demon_Evidence <- function(Comparison) {
  if (!("tidyr" %in% installed.packages())) {
    install.packages("tidyr", dependencies = TRUE)
  }
  if (!("dplyr" %in% installed.packages())) {
    install.packages("dplyr", dependencies = TRUE)
  }

  library(tidyr)
  library(dplyr)

  bands <- unique(Comparison$Band)

  evis <- list()

  for (i in 1:length(bands)) {
    tmp <- Comparison %>%
      filter(Band == bands[i])

    wts <- tmp$AICw
    names(wts) <- tmp$Model

    MxWeight <- max(wts)
    evi <- data.frame(MxWeight / wts)
    evi$Model <- c("BM", "OU", "EB")
    evi$Band <- bands[i]
    evi[evi == 1] <- NA
    names(evi) <- c("ER", "Model", "Band")
    evis[[i]] <- evi
  }

  evis <- do.call(rbind, evis)
  rownames(evis) <- NULL

  evis <- evis[, c(2, 3, 1)]

  return(evis)
}

##### Function to run phylogenetic signal using spectra #####

### Assuming measurement error

demon_K_ME <- function(spectra, spec_ME, tree,
                       nSIM = 1000, nBands = 10) {
  source("https://raw.githubusercontent.com/jesusNPL/PhyloSignal/master/demon_K_test_ME.R")

  ## Inits
  nBands <- nBands
  taxa <- spectra[, 1]
  spectrum <- spectra[, 2:(nBands + 1)]
  spec_ME <- spec_ME[, 2:(nBands + 1)]

  bandNames <- names(spectrum)[1:nBands]

  phylosig_K_ME <- list()

  for (i in 1:nBands) {
    # svMisc::progress(i, max.value = nBands)

    band <- setNames(spectrum[, i], taxa)
    ME <- setNames(spec_ME[, i], taxa)

    tmp_ME <- demon_K_test_ME(
      tree = tree, trait = band, ME = ME, method = "K",
      test = TRUE, nsim = nSIM, bounds_sim = c(-Inf, Inf)
    )

    tmp_ME$Band <- bandNames[i]
    tmp_ME$nSIM <- nSIM

    phylosig_K_ME[[i]] <- tmp_ME

    print(paste0("PS estimated for band ", bandNames[i], " with Measurement Error"))
  }

  phylosig_K_ME <- data.frame(do.call(rbind, phylosig_K_ME))
  rownames(phylosig_K_ME) <- NULL

  phylosig_K_ME <- phylosig_K_ME[, c(10, 11, 1:9)]

  return(phylosig_K_ME)
}

### Assuming NO measurement error
demon_K <- function(spectra, spec_ME, tree,
                    nSIM = 1000, nBands = 10) {
  source("https://raw.githubusercontent.com/jesusNPL/PhyloSignal/master/demon_K_test.R")

  ## Inits
  nBands <- nBands
  taxa <- spectra[, 1]
  spectrum <- spectra[, 2:(nBands + 1)]
  bandNames <- names(spectrum)[1:nBands]

  phylosig_K <- list()

  for (i in 1:nBands) {
    # svMisc::progress(i, max.value = nBands)

    band <- setNames(spectrum[, i], taxa)

    tmp_ME <- demon_K_test(
      tree = tree, trait = band, method = "K",
      test = TRUE, nsim = nSIM, bounds_sim = c(-Inf, Inf)
    )

    tmp_ME$Band <- bandNames[i]
    tmp_ME$nSIM <- nSIM

    phylosig_K[[i]] <- tmp_ME

    print(paste0("PS estimated for band ", bandNames[i], " without Measurement Error"))
  }

  phylosig_K <- data.frame(do.call(rbind, phylosig_K))
  rownames(phylosig_K) <- NULL

  phylosig_K <- phylosig_K[, c(10, 11, 1:9)]

  return(phylosig_K)
}
