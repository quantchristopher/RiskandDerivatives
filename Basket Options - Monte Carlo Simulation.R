##### Pricing Index Options via Naive Monte Carlo Basket Method by Christopher Thomson

### Introduction ####

### This script intends to price a basket or index option on N - number of
##  underlying risky securities. All underlyings are assumed to follow Geometric  
##  Brownian Motion random walks with some degree of correlation between each 
##  respective constituent. Correlation is assumed to be a deterministically static
##  quantity that is derived from the measurable covariance of past time series distributions 
##  between each constituent. From the covariance matrix, the Cholesky Decomposition of the 
##  correlation matrix cross multiplied with random standard normal variables drives is used
##  to generate random walks through successive simulations. A final expectation distribution 
##  is elucidated following 10000-50000 iterations. For purposes of illustration, an example 
##  dataset contained within the qrmdata package, "DJ_Const" consists time series price
##  data for all Dow-Jones Constituent stocks from ~1960-2015. Options will be priced
##  on a basket of certain stocks; in this case a basket of "tech" stocks.

### Data Setup ####
library(qrmdata)
library(qrmtools)

## Import Data:
data("DJ_const")

## Trim to Timeline:
Prices.10s <- DJ_const['2010-01-01/2012-12-31'] # 3 Years worth of pricing data for all constituents

## Create a basket:
Tech.10s <- Prices.10s[,c("AAPL","GE","IBM","INTC","MSFT","VZ")] # 6 Tech Stocks

## Compute Log-Returns:
X.Tech <- returns(Tech.10s, method = "logarithmic")

## Identify the correlation matrix:
COV.Tech <- cov(X.Tech) # Generates a 6x6 covariance matrix from the log-returns:
CORR.Tech <- cov2cor(COV.Tech) # Generates the corresponding correlation matrix

## Isolate cholesky factor of matrix (needed for simulation): 
Chol.Tech <- chol(CORR.Tech)
Lower.Chol.Tech <- t(Chol.Tech) # *** NEED Lower triangle matrix 


## Extract individual variances of each stock:
vars.Tech <- diag(COV.Tech)
sigma.Tech <- sqrt(vars.Tech)

### Generating Price Paths for Pricing a European Call Option ####
### In order to price a European call option on the basket of tech stocks above, need 
#   to generate "price paths" for each constituent stock. Realistically, because this 
#   is a european/vanilla option, there is no need to focus on the actual path and only
#   on the stock's terminal value is of interest. Under the properties of GBM, the 
#   stochastic term is given by Wiener Process W(t). Assuming that there are 6 stocks,
#   there are 6 processes where W_i(t) ~ N(0,t) for i = 1 to 6 and where processes 
#   W_i(t) and W_j(t) have correlation rho_ij. As computed above, let CORR.Tech be the
#   correlation matrix representing this relationship, which can be decomposed into its
#   cholesky factorization Chol.Tech. Because W_1(T),...,W_6(T) follow a multivariate normal
#   distribution ~ N.vec(0.vec,SIGMA) where SIGMA = [diag_6(sqrt(T))]*[CORR.Tech]*diag_6(sqrt(T)),
#   terminal values for each constituent stock price can be simulated via sqrt(T)*Chol.Tech*Norm.Vec.


## Begin by demonstrating a single simulation step:
dim(Tech.10s) # Number will be used to identify current price of each constituent stock from data 
S.T_i <- matrix(NA, nrow = 1, ncol = 6, byrow = FALSE) # Dummy entry matrix to be filled by simulation iteration
Norm.Vec <- rnorm(6, mean = 0, sd = 1) # Naive number generator - randomized

## Simulating the differential dS_i = mu*dt + sigma*dW_i integrated over [0,T] for i = 1,..,6
for(i in 1:6)
  {
    S.T_i[1,i] <-  Tech.10s[754,i]*exp((0.03 - 0.5*sigma.Tech[i]^2)*1 + 
                                         sigma.Tech[i]*sqrt(1)*(Lower.Chol.Tech%*%Norm.Vec)[i]) 
}

B.T <- exp(-0.03*1)*max((sum(S.T_i) - 350),0) # 



### Naive European Basket Call as a Function ####
#   n    = Number of underlying securities
#   S.0s = Initial Prices of each security (must be entered as a vector of n entries)
#   K    = Strike Price of the option
#   r    = Risk Free Interest Rate/Discount Rate
#   COV  = Covariance Matrix of underlying securities (must be entered as nxn matrix)
#   T    = Maturity (Typically in terms of years)
#   NOS  = Number of Simulations 

EuroCall_Basket.MC <- function(n, S.0s, K, r, COV.matrix, T, NOS)
{
  
  Uni.vars    <- diag(COV.matrix) # Extracts the univariate variance of each security's log return
  Uni.sigmas  <- sqrt(Uni.vars) # Converts to standard deviation
  CORR.matrix <- cov2cor(COV.matrix) # Computes the direct correlation matrix 
  Chol.lower  <- t(chol(CORR.matrix)) # Computes the Cholesky Decomposition of CORR (for simulation)
  Normal.vec  <- rnorm(n, mean = 0, sd = 1) # Random 
  
  Low_x_Norm  <- Chol.lower%*%Normal.vec # Key component of simulating final underlying prices; captures correlated Brownian Motion 
  # of each underlying.
  
  Sim.Vec <- rep(NA, NOS) # Dummy Entry Vector
  
  for(m in 1:NOS) # For number of simulations requested (NOS)
    {
      S.T <- rep(NA, n)
      for(i in 1:n) # For number of assets in basket
        {
          S.T[i] <- S.0s[i]*exp((r - 0.5*Uni.sigmas[i]^2)*T + Uni.sigmas[i]*sqrt(T)*Low_x_Norm[i])
        } # End inner for-loop
       Sim.Vec[m] <- exp(-r*T)*max((sum(S.T) - K),0) # Computes the value of the call given strike K
    } # End outter for-loop
     
   1/NOS*sum(Sim.Vec) # Computes average of all basket calls simulated - This is an expectation value that is deemed a sufficient for price for a basket option
}

EuroCall_Basket.MC(n = 6, S.0s = as.vector(c(71.41594,
                                           18.98282,
                                           177.6485,
                                           18.65815,
                                           24.56576,
                                           37.81606)), K = 350, r = 0.03, COV.matrix = COV.Tech, T = 1, NOS = 1000)



