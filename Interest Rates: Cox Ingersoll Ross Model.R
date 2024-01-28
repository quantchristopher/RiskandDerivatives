##### Single Factor Stochastic Interest Rate Simulation - Cox Ingersoll Ross Model

## Following code explores and simulates of the CIR Single Factor Stochastic Interest Rate Model

### Like the Vasicek model, the CIR model assumes a mean reverting Ornstien-Uhlenbeck process.
#   However unlike Vasicek, the stochasticity of the rate is bounded by the interest 
#   rate itself. See the SDE:
#
#                   dr_t = kappa(theta-r_t)dt + sigma*sqrt(r_t)*dW_t
#
#   Note that the risk-neutral drift coffecient is the same as Vasicek, and nearly the same 
#   stochastic scaling coefficient with the exception of sqrt(r_t). 

### Cox Ingersoll Ross Simulation:

# As with Vasicek, can use a discretized process to simulate instantaneous changes to the spot
# rate. The size of the time step between observances will determine the scaling of the 
# parameters in both the drift and stochastic term. For simplicity this simulation assumes the
# impact of the time step has already been incorporated for parameters. However, if this were 
# not implementable, the process would be of the form:
#
#                 r_(t+d.t) = r_t + kappa(theta - r_t)*d.t + sigma*sqrt(r_t*d.t)*dW_t
#
# where d.t is the time step increment.The function below assumes d.t is already incorporated.

# Note that obs is the number of observations of the process.

CIR.Sim <- function(r.0, kappa, theta, sigma, obs){
  CIR.vec <- rep(NA, obs)
  CIR.vec[1] <- r.0
  for(i in 2:obs){
    CIR.vec[i] <- CIR.vec[i-1] + (kappa*(theta - CIR.vec[i-1]) + sigma*sqrt(CIR.vec[i-1])*rnorm(1))
  } # End for-loop
  CIR.vec
} # End Function

## Simulating a single rate curve:
X <- CIR.Sim(r.0 = 0.05, kappa = 0.09, theta = 0.065, sigma = 0.0003, obs = 1260)
matplot(X, type = "l")

## Simulating multiple curves withe same parameters (no-seed):
Trajectories.CIR <- replicate(10, CIR.Sim(r.0 = 0.0666, kappa = 0.07, theta = 0.065, sigma = 0.004, obs = 1500))
matplot(Trajectories.CIR, type = "l", xlab = "Daily Time", ylab = "Observed Spot Rate", 
        main = "1500 Simulated Spot Rates for Given Parameters")

## Note the specifics of each parameter and how it affects the simulation. Kappa determines
#  the rate of mean reversion as t progresses -- large kappa implies "fast" convergence to 
#  a steady state mean. Theta is the value of the long-term steady state mean that the interest
#  is expected to converge to. Sigma regulates the volatility of the "next rate observance". 
#  Somewhat obviously, this means high kappa, low sigma will force a simulated rate to theta 
#  rapidly. Inversely, small kappa, high sigma will exhibit little if any true convergence. 
#  Also note, time scaling is key (something that is not incorporated into the above function). 

## Toggling kappa values:
Trajectories.Kappa <- lapply(c(0.01, 0.03, 0.05, 0.07, 0.09), 
                             FUN = function(kappa){CIR.Sim(r.0 = 0.0666, kappa, theta = 0.04, sigma = 0.00035, obs = 1260)})
matplot(Trajectories.Kappa[[1]], type = "l", xlab = "Time", ylab = "Spot Rate", 
        main = "Simulated Spot Rate Paths for Varied Values of Kappa")
for(i in 2:5){
  lines(Trajectories.Kappa[[i]], type = "l", col = sample(1:5,1))
}
# Increasing kappa values implies faster steady-state convergence, i.e reaches theta value at smaller
# and smaller values of t


## Toggline sigma values:
Trajectories.Sigma2 <- lapply(c(0.0001, 0.0003, 0.0006, 0.0012, 0.0024, 0.0048), 
                             FUN = function(sigma){CIR.Sim(r.0 = 0.047, kappa = 0.075, theta = 0.055, sigma, obs = 1260)})
matplot(Trajectories.Sigma2[[1]], type = "l", ylim = c(0.045, 0.065), xlab = "Time", ylab = "Spot Rate", 
        main = "Simulated Spot Rate Paths for Varied Values of Sigma")
for(i in 2:6){
  lines(Trajectories.Sigma2[[i]], type = "l", col = sample(8:17,1))
}

## Increasing sigma values implies inability to fully converge. Rate will cycle around theta value
#  but will not converge to a "linear" form when sigma is high. 


### Closed Form Bond Pricing Formula ####

## If the interest rate is believed to follow a CIR rate model with well defined parameters,
#  a Zero Coupon Bond's "fair price" can be calculated using the following formula:

CIR.ZCB <- function(Notional, r.t, kappa, theta, sigma, T, t){
  h     = sqrt(kappa^2 + 2*sigma^2)
  A.tau = (2*h*exp((kappa+h)*(T-t)*0.5)/(2*h+(kappa+h)*(exp(h*(T-t))-1)))^((2*kappa*theta)/sigma^2)
  B.tau = (2*(exp(h*(T-t))-1))/(2*h+(kappa+h)*(exp(h*(T-t))-1))
  
  Notional*A.tau*exp(-r.t*B.tau)
}

CIR.ZCB(Notional = 1000, r.t = 0.059, kappa = 0.09, theta = 0.06, sigma = 0.005, T = 7, t = 0)
