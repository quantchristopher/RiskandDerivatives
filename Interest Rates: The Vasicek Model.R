##### Single Factor Interest Rate Model Simulator - Vasicek by Christopher Thomson

## The following code explores the underpinnings of the stochastic interest rate
## first postulated by Oldrich Vasicek. Simulation of the model's differential
## has been converted into an interative process using the Euler-Discretization formula. 

### Primary Vasicek Model ####

# The code below assumes the infinitesimal interest rate aka spot rate follows
# an Ornstein-Uhlenbeck process of the form:
#             
#             dr_t = kappa(theta - r_t)dt + sigma*dW_t
#
# where kappa(theta - r_t) is the risk-neutral measure drift of the process in 
# sequential time, sigma is the standard deviation of distribution of previous 
# spot rate observances and dW_t is the stochastic "disturbance" given by the 
# increment of a Brownian Motion where W_t ~ N(0,t) i.e. mean = 0 and sd = sqrt(t)

## Vasicek Simulating Process:

# For a simpler simulation in which time step considerations are already incorporated into
# the parameters kappa, theta, sigma:
Vasicek.Sim <- function(r.0, kappa, theta, sigma, obs){ 
  # Note that obs represents the number of desired simulated observations
  Vas.Vec <- rep(NA,obs)
  Vas.Vec[1] <- r.0
  for(i in 2:obs){
    Vas.Vec[i] <- Vas.Vec[i-1] + (kappa*(theta - Vas.Vec[i-1]) + (sigma*rnorm(1, mean = 0, sd = 1)))
  } # End For-Loop
  Vas.Vec
} # End Function

## An example of the function to simulate several trajectories of a daily spot rate over
#  next 5 years (with 252 trading/rate quote days each):
Trajectories.1 <- replicate(10, Vasicek.Sim(r.0 = 0.05, kappa = 0.065, theta = 0.065, sigma = 0.0003, obs = 1260))
matplot(Trajectories.1, type = "l", xlab = "Time", ylab = "Spot Rate", xaxt = "no", main = "Vasicek Model Trajectories")


## Second example illustrate the importance of well-understood parameters. Contrasting with the first example above 
#  the parameters input here causes uncharacteristically odd fluctuations in the spot rate. Consequently, all simulated 
#  spot rate paths appear to lack the ability of convergence to a long term steady state rate, something that should be
#  required of an O-U process. Most notably problematic is the instances of spot rate paths
#  becoming negative, a phenomenon that would be rarely observed in actual financial markets. 

Trajectories.2 <- replicate(10, Vasicek.Sim(r.0 = 0.046, kappa = 0.01, theta = 0.001, sigma = 0.003, obs = 1260))
matplot(Trajectories.2, type = "l")

## This can also be illustrated by varying the values of sigma:

Trajectories.Sigma <- lapply(c(0, 0.0002,0.0004,0.0006), FUN = function(sigma){Vasicek.Sim(r.0 = 0.05, kappa = 0.02, theta = 0.065, sigma, obs = 1260)})

matplot(Trajectories.Sigma[[1]], type = "l", ylim = c(0.04,0.077), xlab = "Time", ylab = "Theoretical Spot Rate", main = "Trajectories of Varied Sigma Values")
lines(Trajectories.Sigma[[2]], type = "l", col = "red")
lines(Trajectories.Sigma[[3]], type = "l", col = "blue2")
lines(Trajectories.Sigma[[4]], type = "l", col = "purple2")

# Note how the increasing values in volatility increases the noise around the steady state
# mean of distribution i.e. it becomes increasingly difficult for the process to converge tp
# a long-term seemingly "non-stochastic" steady state rate. 

## And Varying the values of kappa:

Trajectories.kappa <- lapply(c(0.002, 0.02, 0.2), FUN = function(kappa){Vasicek.Sim(r.0 = 0.05, kappa, theta = 0.065, sigma = 0.0003, obs = 1260)})

matplot(Trajectories.kappa[[1]], type = "l", xlab = "Time", ylab = "Theoretical Spot Rate", main = "Trajectories of Varied Kappa Values")
lines(Trajectories.kappa[[2]], type = "l", col = "blue2")
lines(Trajectories.kappa[[3]], type = "l", col = "orange2")
abline(h = 0.065) # See that the steady state conditional expectation is equal to theta. 
# The larger value of kappa forces the process to converge quicker.

# Note that the larger the kappa value, the more rapid the theoretical rate converges to the 
# conditional expected rate, which is given by theta as maturity increases to infinity. 

### Zero Coupon Bond Pricing Function ####

# If the spot interest rate is assumed to follow the Vasicek model, the fair price of a zero 
# coupon bond can then be found using the following function:

# Notional is the cash amount to be redeemed at maturity.
# t is current time.
# T is maturation time.
# kappa is the mean reversion parameter.
# theta is the mean parameter.
# sigma is the volatility parameter.

Vasicek.ZCB <- function(Notional, r.t, t, T, kappa, theta, sigma)
{
  B.tau <- (1 - exp(-kappa*(T-t)))/kappa
  A.tau <- (B.tau - (T-t))*(theta - sigma^2/(2*kappa^2)) - ((sigma^2*B.tau)/4*kappa)
  
  Notional*exp(A.tau - r.t*B.tau)
}

# Example of implementation
Vasicek.ZCB(1000,0.02,8,10,0.02,0.065,0.0003)


