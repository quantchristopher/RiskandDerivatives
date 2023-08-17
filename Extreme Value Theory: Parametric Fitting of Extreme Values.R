##### Parametric Modeling of Extreme Returns Observed in Dow Jones Constituent Stocks by Christopher Thomson

## The following code intends to extract extreme returns (above a given threshold return level) observed 
## in 12 Dow Jones stocks. Those extreme returns are then used to parametrically fit an Extreme Value Distribution
## as well as a Generalized Pareto Distribution, both used as means to model dramatic price swings that may occur
## in the constituent stocks - from which various risk measures (i.e. VaR and ES) can be calculated.


### Load Packages:
library(xts)
library(qrmdata)
library(qrmtools)

#### Data Exploration ####

data("DJ_const")
## Log-Returns of the Dataset
DJ.logrets2 <- diff(log(DJ_const))[-1,]
DJ.logrets <- returns(DJ_const, method = "logarithmic") # Both compute same dataset

## Aggregate by week
DJ.weekly <- apply.weekly(DJ.logrets, FUN = colSums)
## And convert to weekly percentage losses
DJ.weekpercent <- -(exp(DJ.weekly) - 1) * 100
## Crop data to start of 2000s
DJ.weekpercent <- DJ.weekpercent['2000-01-01/2015-12-31']
plot.zoo(DJ.weekpercent[,1:12], type = "h") ### Generates 12 time series plots of each stock's
## weekly percentage loss meaning negative values are negative losses aka gains

## Generate mean excess plots for all 12 stocks
op <- par(mfrow = c(3,4), mar = c(2,2,2,1))
for (i in 1:12){
  data <- DJ.weekpercent[,i]
  mean_excess_plot(data[data > 0], xlab = "", ylab = "", ## plots data > 0 aka omits gains
                   main = names(DJ_const)[i])
}
par(op)

### Focusing on Apple:
APPL <- DJ.weekpercent[,"AAPL"]
mean_excess_plot(data[data > 0]) # plots losses greater than zero
u6 <- 6 ## threshold for losses greater than 6%
abline(v = u6)

## Effect of changing threshold on xi
GPD_shape_plot(APPL) ### Models how increasing the threshold for loss changes the
## estimate for the shape parameter of a GPD distribution modelling the likelihood of
## excess over a given threshold
abline(v = u6)
abline(h = 0.25, lty = "dotted") # assumed that models above this line will produce
## GPD distributions with infinite 4th central moment aka infinite kurtosis. 

#### Fitting GPD model to excess over threshold ####
APPL.exceed <- APPL[APPL > u6] # identifies weekly percent losses that exceed 6%
APPL.excess <- APPL.exceed - u6 # Net of of excess over thresh
APPL.fitGPD <- fit_GPD_MLE(APPL.excess) # Fitting to GPD for shape and scale pars

## QQ Plot using fitted parameters
qqAPPL <- function(p) (
  qGPD(p, shape = APPL.fitGPD$par[1], scale = APPL.fitGPD$par[2])
)
qq_plot(APPL.excess, FUN = qqAPPL)
qq_plot(APPL.exceed, FUN = function(p) (u6 + qqAPPL(p)))

#### Computing semi-parametric risk measure estimators ####
alpha <- c(0.99, 0.995)
VaR <- VaR_GPDtail(alpha, threshold = u6, p.exceed = mean(APPL > u6),
                   shape = APPL.fitGPD$par[1], scale = APPL.fitGPD$par[2])
