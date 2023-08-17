##### Expected Default Frequency Example by Christopher Thomson

## The following code intends to demonstrate Merton's modified method for estimating default frequency on  
## maturing debt securities. Using Radioshack Inc.'s ($RSHCQ) closing stock price data over a successive 
## period, a statistical distribution of stock price returns are estimated. Parameters of this distribution 
## are then input into the famous Black-Scholes formulaic solution. As per the Merton method, from this approach
## an accurate range of values for Radioshack's corresponding debt's Expected Default Frequency and Distance to 
## Default can be estimated given the debt's time to maturity. See comments below for greater detail:

#### Setup and Data Import ####
## Necessary Packages:
library(xts)
library(qrmdata)
library(qrmtools)

## Data Refinement
data("RSHCQ")
RadioShack <- RSHCQ['2010-01-01/2012-03-31']

## Outstanding equity shares by million (approximately):
N.shares <- 100
S.daily <- RadioShack*N.shares # Approximate daily Market Value of total 
# shareholder equity (in millions)
plot.ts(S.daily)

## Create dummy xts object that is indexed by S.daily:
V.daily <- xts(rep(NA,length(S.daily)), time(S.daily))

## Value of one-year debt in millions approximately:
B <- 1042 # aka 1,0420,000,000 total value of short term debt

## Root finding equation:
rooteqn <- function(V,S,t,r,sigmaV,B,T) 
{
  S - Black_Scholes(t,V,r,sigmaV,B,T,"call") 
}
# Computes Value S = S.daily minus the value of a euro-call on underlying V. 
# For use in the for-loop first iteration below.


## Empirical estimate of volatility of total equity:
S.daily.X <- diff(log(S.daily))[-1] # Computes log returns of aggregate 
# market traded equity
Sigma.V <- as.numeric(sd(S.daily.X))*sqrt(250) # Annualized volatility from
# log returns

## First iteration:
for(i in 1:length(S.daily)) 
{
  tmp <- uniroot(rooteqn, interval = c(S.daily[i],10*S.daily[i]), S = S.daily[i], 
                 t = 0, r = 0.03, sigmaV = Sigma.V, B = B, T = 1)
  V.daily[i] <- tmp$root
} # Recall rooteqn has arguments (V,S,t,r,sigmaV,B,T). Uniroot searches for the value
# of the first argument of a function from a given range. Here, V is set to be solved 
# in range (S.daily[i],10*S.daily[10]) and uniroot seeks the value of V that makes the 
# value of a BSM euro-call on V equal so that S.daily - C^{BS} = 0 hence establishing  
# the true value of the firm's total assets. 

Sigma.V.old <- Sigma.V # Reassign calculated value of Sigma.V above
V.daily.X <- diff(log(V.daily))[-1] # Log-returns on calculated daily total asset values 
Sigma.V.new <- as.numeric(sd(V.daily.X))*sqrt(250) # New annualized volatility

Sigma.V.new-Sigma.V.old

## Further iterations:
count <- 1
while ((abs(Sigma.V.new-Sigma.V.old)/Sigma.V.old) > 0.000001) 
{
  count <- count + 1
  for(i in 1:length(S.daily)) 
  { tmp < uniroot(rooteqn, interval = c(S.daily[i],10*S.daily[i]), S = S.daily, 
                  t = 0, r = 0.03, sigmaV = Sigma.V.new, B = B, T = 1)
    V.daily[i] <- tmp$root
  }
  Sigma.V.old <- Sigma.V.new
  V.daily.X <- diff(log(V.daily))[-1]
  Sigma.V.new <- as.numeric(sd(V.daily.X)*sqrt(250))
  
  if ((abs(Sigma.V.new-Sigma.V.old)/Sigma.V.old) < 0.000001){
    stop("Finished")
  }
    
}

Sigma.V.new
count
tail(V.daily)
plot(V.daily, ylim = range(V.daily, S.daily)) 

## Merton Distance to Default:
DD <- (log(V.daily[length(V.daily)]) - log(B)/Sigma.V.new)

(log(V.daily['2011-06-08']) - log(B))/Sigma.V.new

