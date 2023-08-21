##### Duan's GARCH Model for Non-Constant Volatility by Christopher Thomson

## This script reproduces Duan's model for GBM price paths by incorporating non-constant/dynamic volatility through use of 
## General Auto Regressive Conditional Heteroscedasticity - GARCH modeling. In other words, the price differential of an asset
## (in the example Microsoft Share Price) is modeled as dS = mu*dt + sigma*dW, however while sigma is deterministic it is not
## static. Conditioning upon past return data, an ARMA(0,0)+GARCH(1,1) is constructed to model the distribution of future returns.
## The distribution's dynamic sigma parameter provides insight into a more precise term structure model for sigma leading to more
## accurate stochastic term sigma*dW for modeling future price paths. 

### Generating a simulated price path from a Garch model: ####
library(rugarch)
library(qrmtools)

# GARCH(1,1) model specification:
Spec.1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                     distribution.model = "norm",
                     fixed.pars = list(omega = 0.02, alpha1 = 0.15, beta1 = 0.80))

# Generating paths from model above:
m <- 2 # Number of paths
n <- 2000 # Sample Size
(Spec.1.Paths <- ugarchpath(Spec.1, n.sim = n, n.start = 50, m.sim = m)) # Note that n.start is a burn-in parameter

# Extracting values from the simulation:
X <- fitted(Spec.1.Paths) # Simulated return series X_t = mu_t + sqrt(h_t)*Z_t
GARCH.sigmas <- sigma(Spec.1.Paths) # volatilities at each time step aka conditional standard deviations

# Plots of the extracted data:
plot(X[,1], type = "l", ylab = expression(X[1]))
plot(X[,2], type = "l", ylab = expression(X[1]))

plot(GARCH.sigmas[,1], type = "h")
plot(GARCH.sigmas[,2], type = "h")


### Conditioning sigma from log-return data: #### 
## Read-in Microsoft stock price data from September 27, 2000 to September 27, 2001
Stock.MSFT <- read.csv.zoo("MSFT.csv", read = read.csv) # Creates time based object

MSFT <- as.data.frame(Stock.MSFT) 
MSFT <- xts(MSFT, order.by = index(Stock.MSFT)) # Create xts object of the data

X <- returns(MSFT$Close, method = "logarithmic") # Calculate log-returns

## Fitting a ARMA(0,0) GARCH(1,1) model to the MSFT prices:
Spec.2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                     distribution.model = "std")

## Repeating the process for NASDAQ composite:
NASDAQ.3yr <- NASDAQ["2010-01-01/2012-12-31"]
NAS.X <- returns(NASDAQ.3yr, method = "logarithmic")
sd(NAS.X)

NASDAQ.Spec.2Fit <- ugarchfit(spec = Spec.2, data = NAS.X)

plot(NASDAQ.Spec.2Fit) ## Will provide multiple plots available to the user

## Generating a theoretical future price path of MSFT prices that incorporates parameters derived from the GARCH(1,1) regression above:
Sim.Spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                       mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                       distribution.model = "norm",
                       fixed.pars = list(omega = 0.00000365, alpha1 = 0.094 , beta1 = 0.885)) 




Sim.Spec.Paths <- ugarchpath(Sim.Spec, n.sim = 504, n.start = 400, m.sim = 5)

NASDAQ.GSigs <- sigma(Sim.Spec.Paths)

