##### Fitting Multivariate Normal and Multivariate t distributions to financial return data:
library(xts)
library(nvmix) # Package great for fitting distributions (see Ruppert for alternative packages)
library(qrmdata)
library(qrmtools)

#### Data Extraction and Setup ####
data("DJ_const")
str(DJ_const)

## Cutting data down to 5 stock returns in the new millennium:
S <- DJ_const['2000-01-01/', c("AAPL", "BA", "INTC","IBM","NKE")]
## Plotting 5:
plot.zoo(S, xlab = "Time Daily", ylab = c("Apple", "Boeing", "Intel", "IBM", "Nike"), main = "Daily Price Evolution", col = "red4")

## Daily Losses i.e. -log returns
X <- -returns(S, method = "logarithmic") # Computes matrix of -log returns for each stock
pairs(as.matrix(X), main = "Scatter Plot Matrix of Risk-Factor Changes", gap = 0, pch = ".")
# Creates mega scatter plot of covariance/correlation 
plot.zoo(X, xlab = "Daily", main = "")
dim(X)
# Effectively each instance X_t is a 5 component random vector that produces 4024 
# observations which are assumed to be iid in practice. 

#### Fitting and Simulation ####
mu <- colMeans(X) # Average loss for each stock stored in single vector
Sigma <- cov(X) # Covariance matrix of losses
stopifnot(all.equal(Sigma, var(X)))
n <- nrow(X) # number of rows (daily returns) in X
set.seed(271)
X.norm <- rNorm(n, loc = mu, scale  = Sigma) ## Generates 4024 5-dimensional random vectors
## with distribution mu and Sigma of empirical data above. 

## Fitting a multivariate t distribution to X data
fit.t <- fitStudent(X) # Fitting a multivariate t 
X.t <- rStudent(n, df = fit.t$df, loc = fit.t$loc, scale = fit.t$scale) # Generates 4024 5-d RVecs
# with fitted t_distribution 

## Plotting Original Sample, Fitted Normal Sample, and Fitted t Sample:
data.3 <- rbind(original = as.matrix(X), norm = X.norm, t = X.t) # Creates mega matrix
cols <- rep(c("royalblue","black","maroon3"), each = n) # Creates matrix of character strings of colors
pairs(data.3, gap = 0, pch = ".", col = cols)

plot(data.3[,4:5], col = cols, main = "Losses Modeled") #Recall that maroon is the simulated t, it captures
# extreme losses better than the fitted normal distribution. 
