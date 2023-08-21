##### Computing efficient portfolios of N assets using quadratic programming (Not Lagrange)

### Example 16.6 Using CRSPday data
library(Ecdat)
library(quadprog)

data("CRSPday")
rets <- 100*CRSPday[,4:6] ### Creates returns matrix (in percentage rather than decimals)
mean_vect <- apply(rets, 2, mean)
cov_mat <- cov(rets)
sd_vect <- sqrt(diag(cov_mat)) ### Extracts the ii-diagonals from COV(rets), sqrts them to get
## standard deviation of each marginal return distribution

Amat <- cbind(rep(1,3), mean_vect) ### Sets the constraints matrix
muP <- seq(0.05,0.14, length = 300) ### Creates target expected portfolio returns between 5% and 14%

sdP <- muP ### Set up storage for stdev's of portfolio returns given varying Er
wts <- matrix(0, nrow = 300, ncol = 3) # Storage for weights

for (i in 1:length(muP))
{
  bvec = c(1,muP[i]) ### Constraint vector b
  result =
    solve.QP(Dmat = 2*cov_mat, dvec = rep(0,3),
             Amat = Amat, bvec = bvec, meq = 2)
  sdP[i] = sqrt(result$value)
  wts[i,] = result$solution
}

plot(sdP, muP, type = "l", xlim = c(0,2.5), ylim = c(0,0.15), lty = 3) ### Tracing the Efficient frontier
mufree <- (1.3/253)
points(0, mufree, cex = 1, pch = "*") ### Shows risk free asset on frontier plot


sharpe <- (muP - mufree/sdP) ### Computes the Sharpe ratios for each combination of w wts
ind <- (sharpe == max(sharpe)) ### Finds max sharpe ratio

wts[ind, ] # Prints the weight combination that gives largest sharpe ratio i.e. tangency portfolio
### Note that in this example, shortselling two provides the most optimal minimized variance






