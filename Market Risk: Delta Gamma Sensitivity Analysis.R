###### Computing Delta-Gamma and Approximations for Sensitivity Analysis to Adverse Movements by a Portfolio
#  by Christopher Thomson

## The following code intends to construct linear (delta) and quadratic (delta-gamma) approximating functions 
## for changes in a portfolio's aggregate value. Based on the previous work of Paul Embrechts and Marius Hoffert,
## the code follows construction of a linear/quadratic loss operator function given financial return datasets. 
## The delta-gamma loss operator function can be considered as a exploratory form of sensitivity analysis for 
## a portfolio (in this case of equities).

### Data + Setup ####

library(xts)
library(qrmtools)
library(qrmdata)

## 
data("SP500_const")
data("VIX")

SP500.90on <- SP500_const['1990-01-01/'] # Computes new dataset of SP500 from Jan 1, 1990 on
SP5.90on.Aggr <- rowSums(SP500.90on, na.rm = TRUE) # Because only constituent prices of the
# index are listed, needed to compute rowsums. Neatly produces 6553 daily prices

## Log-Returns/Bivariate Dataset
SP5.X <- diff(log(SP5.90on.Aggr))
VIX.X <- returns(VIX/100, method = "diff")['1990-01-01/'] # Since VIX is a percentage value,
# need to use "diff" method
X. <- cbind(SP5.X,VIX.X)
any(is.na(X.))
X <- as.matrix(na.fill(X., fill = "extend"))
colnames(X) <- c("SP500", "VIX")
plot(X)

#### Losses: Linear and Quadratic Approximations ####

## Euro-Call Options:
T <- 1 # 1 year maturity
K <- 100 # Strike
S <- 100 # Current Spot Price
r <- 0.03 # Risk Free Rate 
sigma <- 0.25 # Annualized volatility

## Full revaluation with risk factor changes
daily.t <- 1/250
V.t0 <- Black_Scholes(0, S = S, r = r, sigma = sigma, K = K, T = T, type = "call")
V.t1 <- Black_Scholes(daily.t, S = exp(log(S) + X[,"SP500"]),
                      r = r,
                      sigma = exp(log(sigma) +X[,"VIX"]),
                      K = K, 
                      T = T,
                      type = "call")

L <- -(V.t1- V.t0) # Losses based on the given risk factor change in price of underlying
# and risk factor change in volatility from empirical values of S&P500 log returns and 
# log returns on VIX

hist(L, breaks = 100, col = "red", probability = TRUE); box()

## Delta Approximation (linear loss operator)
(Greeks <- Black_Scholes_Greeks(0, S = S, r = r, sigma = sigma, K = K, T = T))
# Recall that the linear loss operator relies on partial derivatives which are
# computed as the BSM Greeks. Therefore, the linear loss operator is given by:
L.linears <- -(Greeks[["delta"]] * X[,"SP500"]*S +
                 Greeks[["theta"]] * daily.t +
                 Greeks[["vega"]] * X[,"VIX"] * sigma)
plot(L, L.linears)
abline(0, 1, col = "royalblue1")


#### A More Explicit Example of Loss Operator on Euro-Options ####
## (Recommended to clear data from above before running this section of script
t <- 0 # Current time
T <- 1 # Time to maturity (in years for BSM function)
K <- 100 # Strike Price
S <- 110 # Stock price now 
r <- 0.02 # Interest rate 2%  
sigma <- 0.2 # Annualized Volatility of returns (assuming calculated from some 
# parametric distribution and/or GARCH)

## Applied for loss operators
Change.S <- 0.05 # change in stock price (5% daily log return)
Change.sigma <- 0.02 # change in implied volatility 
Change.t <- 1/250 # If T = 250/250 
Change.rho <- 0.001 
(Greeks <- Black_Scholes_Greeks(0, S = S, r = r, sigma = sigma, K = K, T = T))

## At t = 0, S_0 create a hedged portfolio where one is short the call option, 
# long the stock at BS-Delta Number of shares. Note that only the option's delta
# has only been hedged. A full replicating portfolio would also require money market
# trades. 
V.t0 <- S * Greeks[["delta"]] - Black_Scholes(t, S = S, r = r, sigma = sigma, K = K, T = T)
# No money market account
S * Greeks[["delta"]] # Value of delta hedge given initial values
Black_Scholes(t, S = S, r = r, sigma = sigma, K = K, T = T) # BSM value of vanilla option
# Risk factor change in the stock price and in the implied annualized volatility
S.new <- exp(log(S) + Change.S)
V.t1 <- S.new * Greeks[["delta"]] - Black_Scholes(t + Change.t, S = S.new, r = r,
                                                  sigma = sigma + Change.sigma, 
                                                  K = K, T = T)
# Loss on portfolio given by:
L1 <- -(V.t1 - V.t0) # This is the implicit change in net value of the long shares
# short call option when the underlying's price changes and the implied volatility changes

# Computing the linear loss operator i.e. l_delta[t](x) aka delta approximation
# ldelta_[t](x) = C_tau * change.t + C_sigma * Change.sigma # Recall that delta term
# goes to zero:
(L.linear <- Greeks[["theta"]]*Change.t + Greeks[["vega"]]*Change.sigma)


