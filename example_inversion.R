
##### Inversion spectra #####

require(PEcAnRTM)

### Prepare data
#source("https://raw.githubusercontent.com/jesusNPL/SpecEvolution/main/R/readSPEC.R")
#state <- "OH"
#ruta <- "/Users/jesusnpl/Dropbox/Oak.Project.SVC data/*Ohio/" 
#metaDT <- read.csv("/Users/jesusnpl/Dropbox/Oak.Project.SVC data/2019.META.DATA/2019.SVC.DATA.LABELS.csv")

#spec_OH <- demon_readSPEC(path = ruta, 
 #                         metadata = metaDT, 
  #                        state = state, 
   #                       format = "sig")

#tt <- spec_OH[[2]]
#spec_splice <- guess_splice_at(spec_OH[[1]])
#[1]  984.0667 1896.3200

#spec_matched <- match_sensors(x = spec_OH[[1]], splice_at = spec_splice,
                                 #interpolate_wvl = c(5, 1))

#spec <- spectrolab::resample(spec_matched, new_bands = seq(400, 2500, 1))
#dta4 <- as.data.frame(spec, fix_names = "none", metadata = FALSE)

load("~/Documents/GitHub/SpecEvolution/Data/test_alba.RData")

wls <- 400:2500

waves <- paste0("Wave_", wls)

tail(waves)

alba <- as.numeric(dta4[1, ][2:2102])
names(alba) <- waves

plot(wls, alba, type = "l", col = "red", lwd = 2, 
     xlab = "Wavelength (nm)", ylab = "Reflectance")

### Params inversion

params <- c(
  "N" = 1.4, "Cab" = 40, "Car" = 15,
  "Cbrown" = 0.5, "Cw" = 0.002, "Cm" = 0.004
)

invert.options <- default.settings.prospect
invert.options$model

invert.options$n.tries <- 4      # Number of attempts
invert.options$nchains <- 2      # Number of MCMC chains
invert.options$ngibbs <- 10000    # Number of iterations per chain
invert.options$burnin <- 2000    # Length of burnin period
invert.options$do.lsq.first <- TRUE # Initialize with results from a fast least-squares optimization algorithm

names(invert.options)

### Run inversion

inversion.alba <- PEcAnRTM::invert.auto(observed = alba,
                                invert.options = invert.options,
                                quiet = TRUE)

### Check convergence
par(mfrow = c(2, 1))
plot(inversion.alba$samples, auto.layout=FALSE)
par(mfrow = c(1, 1))

samples.mat.alba <- as.matrix(inversion.alba$samples)[-(2000:0),1:5]
colnames(samples.mat.alba) <- params.prospect5
pairs(samples.mat.alba, pch = ".")

### Extract traits
means.alba <- unlist(inversion.alba$results[grep("mu", 
                                              names(inversion.alba$results))])[1:5]

means.alba
#N.mu       Cab.mu       Car.mu        Cw.mu        Cm.mu 
#1.352356757 49.973335947  8.514526977  0.005992937  0.005133357 

median.alba <- unlist(inversion.alba$results[grep("med", 
                                                 names(inversion.alba$results))])[1:5]

median.alba
#N.med      Cab.med      Car.med       Cw.med       Cm.med 
#1.354597609 49.635895688  8.171001474  0.005977718  0.005195701 

### Simulate spectra using traits from inversion

prospect.sim.alba <- prospect(means.alba, 5)[, 1]  # reflectance
prospect.sim.median.alba <- prospect(median.alba, 5)[, 1]  # reflectance

#### Plot observed and simulated spectra
plot(wls, alba, type = 'l', col = 1, lwd = 3,  
     xlab = "wavelength (nm)", ylab = "reflectance")
lines(wls, prospect.sim.alba, type = 'l', col = 2, lwd = 3)
lines(wls, prospect.sim.median.alba, type = 'l', col = 3, lwd = 1)
legend("topright", c("observed", "predicted1", "predicted2"), lty = 1, col = 1:3)

save(means.alba, median.alba, inversion.alba, 
     file = "Data/inversion_alba.rda")
