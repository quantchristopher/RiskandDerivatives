##### QRM Chapter 5 - Danish Fire Data

##### Packages + Setup #####

library(xts)
library(qrmdata)
library(qrmtools)

## Danish Fire Data
data("fire")
plot.zoo(fire, xlab = "Year", ylab = "Insurance Claim Losses per 1M DKK", main = "Danish Fire Insurance Claims")

## Sample Mean Excess Plot (with two plausible thresholds)
mean_excess_plot(fire) ## Plots data above a minimal threshold in which the mean excess plot is roughly linear
u10 <- 10 # threshold 1
u20 <- 20 # threshold 2
abline(v = c(u10,u20))

GPD_shape_plot(fire)
abline(v = c(u10,u20))
abline(h = 0.5, lty = 3) ## Models above this line will generally have infinite variance

#### Fitting GPD parameters for a given threshold ####
## For u = 10
exceed_u10 <- fire[fire > u10] # Loss of over 10 million DKK
excess_u10 <- exceed_u10 - u10 # Net excess over threshold
fitted_u10 <- fit_GPD_MLE(excess_u10) # Computes parameters of GPD(xi,beta) aka shape and scale
# via MLE log-likelihood
GPD_u10_par <- fitted_u10$par # shape and scale
print(GPD_u10_par)

## For u = 20 million DKK (in units of 1 million DKK)
exceed_u20 <- fire[fire > u20]
excess_u20 <- exceed_u20 - u20
fitted_u20 <- fit_GPD_MLE(excess_u20)
GPD_u20_par <- fitted_u20$par

#### Plot Comparison of empirical excess distribution vs fitted GPD G_{xi,beta}  ####
## u = 10
empiric_ex <- edf_plot(excess_u10) # Plots empirical distribution of excess 
# stepwise over threshold 10
tail_empiric_excess <- tail(sort(as.numeric(excess_u10)), n = -1)
lines(tail_empiric_excess, pGPD(tail_empiric_excess, shape = GPD_u10_par[1], scale = GPD_u10_par[2]))



#### Semi-Parametric Risk Measure Estimators ####
## VaR_alpha and ES_alpha for two CI alphas (99% and 99.5%) and both thresholds (10 and 20 million DKK)
alpha <- c(0.99, 0.995)
VaR_u10 <- VaR_GPDtail(alpha, threshold = u10, p.exceed = mean(fire > u10),
                       shape = GPD_u10_par[1], scale = GPD_u10_par[2])


