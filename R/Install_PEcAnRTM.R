# These packages are not necessary
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/data.land") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/data.atmosphere") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "models/ed") # check

# In order to install PEcAnRTM in your computer you should install the next 
# packages in this specific order... 

# Also, PEcAnRTM only works in unix like OS, so, don't waste your time 
# in trying to install it in windows. I know, it's a pain, but it is what it is

devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "base/logger") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "base/remote") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "base/utils") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "base/db") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "base/settings") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/benchmark") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/emulator") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/meta.analysis") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/priors") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/uncertainty") # check
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/assim.batch") # check
# Finally
devtools::install_github("pecanproject/pecan", ref = "develop", 
                         subdir = "modules/rtm") # check

## Example 
library(PEcAnRTM)

wl <- 400:2500
params <- c(
  "N" = 1.4, "Cab" = 40, "Car" = 15,
  "Cbrown" = 0.5, "Cw" = 0.002, "Cm" = 0.004
)

p4 <- prospect(params[c(-3, -4)], version = 4)
p4
class(p4)
p4[1:50, 1]
p4[[500:520, 2]]

p5 <- prospect(params[-4], version = 5)
p5b <- prospect(params, version = "5B")
p_multi <- cbind(p4, p5, p5b)
matplot(p_multi[, c(1, 3, 5)], lty = 1:3, col = 2, ylim = c(0, 1))
matplot(1 - p_multi[, c(2, 4, 6)], lty = 1:3, col = 3,  add = TRUE)
legend("topright", c("Reflectance", "Transmittance"), col=c(2, 3), lty = 1)
legend("top", c("4", "5", "5B"), lty = 1:3)

sail.params <- defparam("pro4sail")
print(sail.params)
p4s <- pro4sail(sail.params)
matplot(p4s, xlab = "Wavelength (nm)", ylab = "Reflectance")
legend("topright", as.character(1:4), col = 1:4, lty = 1:4)

