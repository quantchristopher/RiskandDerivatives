### Ch_5 Fitting Distribution by Max-Likelihood
### Example finds a distribution of changes in the risk-free rate 

rm(list = ls())
data(Capm, package = "Ecdat") ### Loads data stored in the Ecdat package, specifically Capm dataset
rf <- Capm$rf
rf_diff <- diff(rf)

### Requires MASS
library(MASS)

fitdistr(rf_diff, densfun = "t") ### Fits data to distribution, in this case t
### Estimated values of t-distribution of rf data are mu = 0.00122, sigma = 0.0458, df = 3.3367
### Numbers in parentheses below are SEs of each value

### Comparison section relies on fGarch (specifically the optim() function for optimizing log-likelihood)
library(fGarch)

n <- length(rf_diff)
start <- c(mean(rf_diff), sd(rf_diff), 5)

dumbfun <- function(x) x+2 ### Reminder of how to write a simple function in R
dumbfun(9) ### Prints 11

### This code returns the negative log-likelihood of the t-density function with given rf data
loglike_t <- function(beta) sum( -dt((rf_diff-beta[1])/beta[2],
                                     beta[3], log = TRUE) + log(beta[2]) )                        


### Optimizes the values of the density function                                 
fit_t <- optim(start, loglike_t, hessian = TRUE, method = "L-BFGS-B",
               lower = c(-1,.001,1) )
### Note that it returns very similar values of mu, sigma, and df as the fitdistribution function

AIC_t <- 2*fit_t$value+2*3
BIC_t <- 2*fit_t$value+log(n)*3
sd_t <- sqrt(diag(solve(fit_t$hessian)))


### The above can be redone more simply using the standardized t-distribution
loglik_std_t <- function(beta) sum(-dstd(rf_diff, mean = beta[1],
                                         sd = beta[2], nu = beta[3], log = TRUE))

fit_std_t <- optim(start, loglik_std_t, hessian = TRUE, method = "L-BFGS-B",
                   lower = c(-0.1, .001, 2.1))
print(fit_std_t)


#### Estimating MLE of Ferrari Return Data
Race <- read.csv("RACE.csv")
logret_ferr <- diff(log(Race$Adj.Close))

### Use density() to estimate KDE 
density(logret_ferr, adjust = 0.5)
plot(density(logret_ferr, adjust = 1.5))

ferr_mu <- mean(logret_ferr)
ferr_sigma <- sqrt(var(logret_ferr))

fitdistr(logret_ferr, densfun = "normal")

### Attempting to write a likelihood function based on the above example. Density is normal here 
loglike_ferr <- function(ferr_mu,ferr_sigma) sum( -log((1/ferr_sigma*sqrt(2*pi))*exp(0.5*((logret_ferr-ferr_mu)/ferr_sigma)^2
                                                                   )) )
loglike_ferr2 <- function(ferr_mu, ferr_sigma) sum(-log(dnorm(logret_ferr, mean = ferr_mu, sd = ferr_sigma )))

### Not finished
start2 <- c(ferr_mu,ferr_sigma,5)
fit_ferr <- optim(start2, loglike_ferr, hessian = TRUE, method = "L-BFGS-B", lower = c(-1, .001, 1))




