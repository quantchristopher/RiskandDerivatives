#### Applying Block Maxima Method to S&P 500 Returns by Christopher Thomson

## The following code is an example of applying Extreme Value Theory, specifically the Block Maxima Method 
## to identify an equity index's extreme return level and extreme return period. Here Standard and Poor's 
## Blue Chip 500 Index equity data is used. 

##  ** Note log-returns are used **

### Setup ####
library(xts) # for functions around time series objects
library(qrmdata) # for the S&P 500 data
library(qrmtools) # for returns()

### Data Importing and Preliminary Analysis ####
data("SP500")
S <- SP500
X <- returns(S, method = "logarithmic") # Computes Log-Returns
L <- -X # Losses in standard loss notation
stopifnot(all.equal(L, -diff(log(S))[-1], check.attributes = FALSE)) ### No idea what this does

## Looking at individual dates
# Friday October 16th 1987
L['1987-10-16'] # Loss was ~0.0529 i.e. SP log returned -5.29%
# Monday October 19th 1987 aka Black Monday
L['1987-10-19'] # Loss was ~0.22899 i.e. SP log returned -22.899% in value

## Consider classical returns i.e. Y_t = (S_t/S_t-1) - 1
## Uniquely, S_t = (1 + beta) * S_{t-1} => Y_t = -(S_t/S_{t-1}-1) = -beta
## => The negative classical returns Y_t give exactly the drop beta (= -beta)
Y <- -returns(S, method = "simple") # classical negative returns
stopifnot(all.equal(Y, -diff(S)[-1]/as.numeric(S[-length(S)]),
                    check.attributes = FALSE))

Y['1987-10-16'] # ~=  5.16% (drop)
Y['1987-10-19'] # ~= 20.47% (drop)


##### Applying Block Maxima Method #####
plot.zoo(L, main = "SP500 Negative Log Returns")

### Extracting Maxima
M.years <- period.apply(L, INDEX = endpoints(L, "years"), FUN = max) # Computes yearly Maxima from 1960-1987
#Generating indexing endpoints for each block of returns
endpts <- endpoints(L, "quarters") # end indices for quarters
endpts <- endpts[seq(1, length(endpts), by = 2)] # Converting endpoints to half years, as data is by quarters
M.halfyears <- period.apply(L, INDEX = endpts, FUN = max) # Computes half-yearly maxima i.e. largest loss (most negative log return risk factor)

## Extracting Maxima for returns between January 1, 1960 and October 16, 1987
Lcrop <- L['1960-01-01/1987-10-16']
M.Lcropyears <- period.apply(Lcrop, INDEX = endpoints(Lcrop,"years"), FUN = max)
endptscrop <- endpoints(Lcrop, "quarters") # Dnd indicies for quarters in crop set
endptscrop <- endptscrop[seq(1, length(endptscrop), by = 2)] # Ending indices for half years in crop set
M.Lcrophalfyears <- period.apply(Lcrop, INDEX = endptscrop, FUN = max)


### Fitting GEV distribution H_{xi, mu, sigma} to both maxima
## Raw Yearly:
GEV.years <- fit_GEV_MLE(M.years)
stopifnot(GEV.years$convergence == 0) # Assures function will stop if no convergence
print(GEV.years)
Xi.years <- GEV.years$par[["shape"]]
Mu.years <- GEV.years$par[["loc"]]
Sig.years <- GEV.years$par[["scale"]]
  
## Raw Half-Yearly:
GEV.Halfyears <- fit_GEV_MLE(M.halfyears)  
stopifnot(GEV.Halfyears$convergence == 0)
print(GEV.Halfyears)
Xi.halfyears <- GEV.Halfyears$par[["shape"]]
Mu.halfyears <- GEV.Halfyears$par[["loc"]]
Sig.halfyears <- GEV.Halfyears$par[["scale"]]  


GEV.cropyears <- fit_GEV_MLE(M.Lcropyears) # Xi ~ 0.2971 > 0 so Frechet Domain
GEV.crophalfyears <- fit_GEV_MLE(M.Lcrophalfyears) # Xi ~ 0.3401 > 0 also Frechet Domain 

##### Computing Exceedence Probabilities #####
## Probability that next year's maximal risk factor change (i.e. largest loss) exceeds
## all previous ones in the cropped data:

GEV.cropyear_pars <- GEV.cropyears$par # Extracts the parameter values of the fitted GEV distribution
## pGEV computes probability of GEV distribution of given quantile with mapped parameters
1 - pGEV(max(head(M.Lcropyears), n = - 1), shape = GEV.cropyear_pars[1], loc = GEV.cropyear_pars[2], scale = GEV.cropyear_pars[3])
## 2.58% chance that next years maximal loss exceeds 22% (i.e. log-return is < -.22)
## ***Note that pGEV() computes the CDF of a given quantile where losses are measured positively
##    This is interpred: pGEV(.10|theta_hat) is cumulative probability the loss is equal to or less than 10% log loss i.e.
##    a log return of -10%, assuming parameters for GEV have been estimated. To compute probability thresholds above, 
##    1-pGEV(.10) is the probability that the loss will exceed 10% i.e. a log-return = -.10. 


## For half year blocks:
GEV.crophalfyear_par <- GEV.crophalfyears$par
1- pGEV(max(head(M.Lcrophalfyears, n = - 1)), shape = GEV.crophalfyear_par[1], loc = GEV.crophalfyear_par[2], scale = GEV.crophalfyear_par[3])
# 1.49% probability the next half year's loss will exceed maximum loss of previous data


### Computing Return Levels for 10 and 50 years i.e. block has n = 260 days
## Recall: k n-block return level = r_{n,k} = H^-1 of (1-1/k) i.e. the (1-1/k)th quantile of H~(xi,mu,sig)
## This level is the value of a loss than can be expected to be exceeded every in out of every
## k n-blocks (so here k = 10, 50 for n = 260[year] or 130[half-year])

## 10 and 50 Year Blocks:
retlevel_10years <- qGEV((1-1/10), shape = GEV.cropyear_pars[1], loc = GEV.cropyear_pars[2], GEV.cropyear_pars[3]) * 100
retlevel_50years <- qGEV((1-1/50), shape = GEV.cropyear_pars[1], loc = GEV.cropyear_pars[2], GEV.cropyear_pars[3]) * 100
## In other words a loss is expected to exceed 4.42% at least once in the coming 10 years
## and a loss is expected to exceed 7.49% at least once in the coming 50 years assuming the 
## parameterized distribution of GEV.cropyears is accurate in explaining maxima of each yearly
## block.

## 20 and 100 Half-Year Blocks:
retlevel_20half <- qGEV((1-1/20), shape = GEV.crophalfyear_par[1], loc = GEV.crophalfyear_par[2], scale = GEV.crophalfyear_par[3])*100
retlevel_100half <- qGEV((1-1/100), shape = GEV.crophalfyear_par[1], loc = GEV.crophalfyear_par[2], scale = GEV.crophalfyear_par[3])*100
## In next 20 semesters, a return is expected to exceed 4.56% and in next 100 semesters a return
## is expected to exceed 7.90% assuming GEV.crophalfyear captures distribution of maxima every
## semester. 
## *** Note for above From the book r_n,k = mu_hat + (sig_hat/xi_hat)((-ln(1-1/k))^-xi_hat - 1)


### Computing return period:
## If H is the density of true distribution of the n-block maximum, the return period
## of the event {M_n > u} is given by k_n,u = 1/barH(u). In other words, in k_n,u is number of 
## blocks we would expect to observe a single block in which the level u is exceeded. If there
## n-sized was a strong tendency for the extreme values to cluster, we might expect to see multiple
## exceedances of the level within the same block - however BMM may not capture these - it will
## only speculate the length of time for a given level, or a given level over a specified 
## length of time. Here, what is the expected return period, based on the calculated GEV,
## of a risk-factor change (i.e. negative return) as large as the one the S&P500 saw on 
## Black Monday? In other words number of n-blocks in which one of them has a return that
## exceeds -22% or loss of 22%.
badmonday <- as.numeric(L['1987-10-19'])

1/(1-pGEV(badmonday, shape = GEV.cropyear_pars[1], loc = GEV.cropyear_pars[2], scale = GEV.cropyear_pars[3]))
# Equals 1875.99 => 1876 future 260 length blocks in which one of them contains a loss of 22%

1/(1-pGEV(badmonday, shape = GEV.crophalfyear_par[1], loc = GEV.crophalfyear_par[2], scale = GEV.crophalfyear_par[3]))
# Equals 2299.93 => 2300 future 130 length blocks in which one of them contains a loss of 22%

### ***Note: that all of these calculations are based on MLE of the GEV distributions parameters.
##     For estimating GEV ~ H(xi, mu, sigma) using profile likelihood, see other R-Script

  
