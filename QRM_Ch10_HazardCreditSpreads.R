##### QRM Chapter 10 - Credit Spreads from Hazard Rates

## Define inputs for both recovery of treasury and recovery of face value models:
delta.LGD <- 0.5 # Calculated loss given default percentage
gamma.con <- 0.04 # Hazard rate function gamma(t), in this case just a constant
epsilon <- 0.01 
rf <- 0.02 # risk free rate function r(t), in this case just a constant
T <- 5 # 5 year maturity


## Recovery of Treasury/Default Free Credit Spread:
RTspread <- function(t,T,gamma.con,delta.LGD)
{
  -log(1-delta.LGD*(1-exp(-gamma.con*(T-t))))/(T-t)
}

t <- seq(from = 0, to = (T-epsilon), length = 50)

C.RT <- RTspread(t = t, T = T, gamma.con = gamma.con, delta.LGD = delta.LGD) 

plot(t, C.RT)

