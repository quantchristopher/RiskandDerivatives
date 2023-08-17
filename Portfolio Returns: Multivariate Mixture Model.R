##### Financial Returns Distribution Modeling by Christopher Thomson

## The following code intends to model a multivariate distribution of returns for a portfolio of Dow Jones
## listed public equities. Using Maximimum Likelihood Estimation, a probability distribution of the mixture
## model variety is elucidated. The distribution can then be used to make inferences about future returns as
## well as provide insights about the portfolio's overall risk profile.

# Requires R Package:
library(ghyp)

#### Data Setup and Extraction ####
data("DJ_const")
d <- 10 
S <- DJ_const['2000-01-01/',1:d] # Selects daily prices of first 10 constituents
# of DJ Index data provided for 2000 on 

## Daily log returns
X.d <- returns(S, method = "logarithmic")
## Weekly log returns
X.w <- apply.weekly(X.d, FUN = colSums)
## Monthly log returns
X.m <- apply.monthly(X.d, FUN = colSums)

#### Uni/Multivariate Tests of Normality ####

### Daily log-returns:
## Shapiro Test:
apply(X.d, MARGIN = 2, FUN = function(x) shapiro.test(x)$p.value)
# Applies the shapiro test to each marginal distribution - specifically returning 
# the p-value. Since all p-values are very small for every constituent, its fair
# to say every stock's return is not normally distributed as extreme values at both
# ends appear present

## Anderson-Darling Test:
maha2_test(X.d) # Using Mahalanobis distances and angles, identifies deviation from 
# normality of joint distribution. Small p-value indicates reject null hypothesis
# indicating that yes, the joint distribution of marginal is not normally distributed

mardia_test(X.d, type = "kurtosis") # Tests for kurtosis of MV data
mardia_test(X.d, type = "skewness") # Tests for skewness of MV data

### Weekly log-returns:
apply(X.w, 2, function(x) shapiro.test(x)$p.value) # All p-values << 0.05, so
# all univariate margins likely non-normally distributed
maha2_test(X.w) 
mardia_test(X.w, type = "kurtosis")
mardia_test(X.w, type = "skewness")

# Literally every p-value calculated is <<< .05 meaning the data is definitely not
# normally distributed 

### Visual Tests
pairs(as.matrix(X.d), gap = 0, pch = ".") # Creates pairwise scatterplots of each
qqnorm(X.d[,1]) # QQNorm plot
qq_plot(X.d[,1], FUN = qnorm, method = "empirical") # QRMTools qq_plot approach

D2.d <- mahalanobis(X.d, center = colMeans(X.d), cov = cov(X.d)) # Computes Mahalonobis Distances
qq_plot(D2.d, FUN = function(p) qchisq(p, df = d)) # Departs from theoretical trendline

### Consider similar visual test for weekly returns as well

#### Fitting Various Mixture Models beyond MVnormal ####
max.iter <- 1e4 # Maximum number of iterations for fitting procedure so that the
# computer is not overloaded 

### Fitting Daily-Log Returns
## Symmetrical Distributions (no skew):
fit.t.sym.d <- fit.tmv(X.d, symmetric = TRUE, nit = max.iter, silent = TRUE) # t via MLE
fit.NIG.sym.d <- fit.NIGmv(X.d, symmetric = TRUE, nit = max.iter, silent = TRUE) # NIG via MLE
fit.GH.sym.d <- fit.ghypmv(X.d, symmetric = TRUE, nit = max.iter, silent = TRUE) # GH via MLE
## Skewed Distributions
fit.t.skw.d <- fit.tmv(X.d, symmetric = FALSE, nit = max.iter, silent = TRUE) # skewed t
fit.NIG.skw.d <- fit.NIGmv(X.d, symmetric = FALSE, nit = max.iter, silent = TRUE) # skewed NIG
fit.GH.skw.d <- fit.ghypmv(X.d, symmetric = FALSE, nit = max.iter, silent = TRUE) # skewed GH

## Comparing log-likelihoods of fit
likelihoods.d <- c(fit.t.sym.d@llh, fit.NIG.sym.d@llh, fit.GH.sym.d@llh,
                   fit.t.skw.d@llh, fit.NIG.skw.d@llh, fit.GH.skw.d@llh)
which.max(likelihoods.d) # 6th value indicates skewed GH has highest llh
