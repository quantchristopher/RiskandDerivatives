##### QRM Chapter 10 - Calibrating a CDS for fixed hazard rate

## General Inputs into CDS:
delta.t <- 0.25 # Quarterly premium payments by protection buyer (t in years)
N <- 20 # Total number of payments for a contract
r <- 0.02 # Risk free rate
x.star <- 0.0042 # Swap Spread 
delta.LGD <- 0.6 # Loss Given Default percentage where notional equals 1

## CDS Leg Calculations:
premium.pmts <- function(gamma,x.star,N,delta.t,r)
  {
  k = 1:N
  x.star*sum(exp(-r*k*delta.t)*delta.t*exp(-gamma*k*delta.t))
} # Calculates the total discounted premium payments

premium.pmts(.04, .0042, 20, 0.25, .02) # Using above inputs and assuming a 4% 
# hazard rate total discounted cash flows made by protection buyer is 0.018 or 
# 1.8% value of the reference entity.

gamma <- seq(from = 0, to = 0.02, length  = 20) # 20 fixed values of hazard rate
pleg <- rep(NA, length(gamma)) # Premium leg dummy vector to be filled with values
# from the below for-loop

for(i in 1:length(gamma))
  {
  pleg[i] <- premium.pmts(gamma[i],x.star, N, delta.t, r)
} # Runs various inputs of gamma into the pricing equation


default.pmt <- function(gamma,delta.LGD,N,delta.t,r)
  {
  delta.LGD*gamma*exp(-(r+gamma)*N*delta.t)/(r+gamma)
} # Function calculates the default leg of the CDS

default.pmt(.02, .6, 20, .25, .02) # Function test


dleg <- rep(NA, length(gamma)) # dummy default leg vector
for(i in 1:length(gamma))
  {
  dleg[i] <- default.pmt(gamma[i], delta.LGD, N, delta.t, r)
} # Fills dummy d-leg vector with the values in of gamma vector


plot(gamma, pleg, ylim = range(dleg,pleg), type = "n")
lines(gamma, dleg, col = "blue3")
lines(gamma, pleg, col = "red2")


root.function <- function(gamma,x.star,N,delta.t,r,delta.LGD)
{
  premium.pmts(gamma,x.star,N,delta.t,r) - default.pmt(gamma,delta.LGD,N,delta.t,r)
}
tail(gamma)

tmp <- uniroot(root.function, interval = c(0,0.02), x.star = x.star, 
               N = N, delta.t = delta.t, r = r, delta.LGD = delta.LGD) 
abline(v = tmp$root)
