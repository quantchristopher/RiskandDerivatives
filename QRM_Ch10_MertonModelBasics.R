##### QRM Chapter 10 - Simulating Basics of Merton Model 

#### Setup and Path Simulation ####

library(qrmtools)
library(Sim.DiffProc)

Brownie1 <- BM(N = 10, M = 1, x0 = 0, t0 = 0) 
Brownie2 <- BM(N = 100, M = 2, x0 = 0, t0 = 0)
Brownie3 <- GBM(N = 364, M = 1, x0 = 1, t0 = 0, T = 1, theta = 0.02, sigma = 0.25)
plot.ts(Brownie3)

## Parameterizing Merton Model:
V.0 <- 1 # Initial Value
mu.V <- 0.03 # Average drift under real world measure
sigma.V <- 0.25 # Annualized measurable volatility
r <- 0.02 # Risk free interest rate
T <- 1 # time to maturity in years
N <- 364 # Observable days
B <- 0.85 # Notional Value of the debt


## Simulate asset value trajectories for Merton model
npaths <- 50
paths <- matrix(NA, nrow = N+1, ncol = npaths)
for (i in 1:npaths) {
  paths[,i] <- GBM(N = N, M = 1, x0 = V.0, t0 = 0, T = T, 
                   theta = (mu.V - r)/sigma.V, sigma = sigma.V)
}


opar <- par(mar = c(3,3,2,1), mgp = c(2,1,0))
times <- (1:(N+1))/(N+1)
plot(times,paths[,1], type = "l", main = "2023 Structural Volatility: Distressed Floor versus Unlimited Cap (Simulated N = 150)",
     xlab = "t", ylab = expression(V[t]), ylim = range(paths))
for(i in 1:npaths) {
  lines(times, paths[,i], col = sample(1:20,1))
}
abline(v = 1)
abline(h = B)

## Overlay a lognormal density:
lnV.mu <- log(V.0) + ((mu.V-r)/sigma.V -0.5*sigma.V^2)*T 
lnV.sigma <- sigma.V*sqrt(T) 
rvals <- range(paths) # Min to max value of all paths
vals <- seq(from = rvals[1], to = rvals[2], length = 100) # 100 unit sequence from min to max of paths
dens <- dlnorm(vals, meanlog = lnV.mu, sdlog = lnV.sigma) # Generates density curve of lognormal
# distribution with quantiles vals and mean/sd calculated from above

dvals <- (max(dens)-dens)/(max(dens)-min(dens))
lines(dvals,vals)


## Below is a path where at T (t = 365th entry), value of V.T < B = 0.85
set.seed(491)
V.t <- GBM(N = N, M = 1, x0 = V.0, t0 = 0, T = T, 
           theta = (mu.V - r)/sigma.V, sigma = sigma.V)
V.t[365] # 0.721 < 0.85. This path would infer a default under the Merton model
V.t[365] - 0.85 # Bondholders would theoretically take control of remaining assets
# shareholder equity is reduced to 0.  
times491 <- seq(from = 0, to = 1, length = N+1)
plot(times491, V.t, type = "l", ylim = range(0.6, max(V.t)), 
     xlab = "Timestep t", ylab = expression(V[t]))
abline(v = 1)
abline(h = 0.85)
# Add default free debt value
lines(times491, B*exp(-r*(T-times491)), col = 3)
# Add defaultable debt
B.t <- B*exp(-r*(T-times491)) - Black_Scholes(times491,V.t,r,sigma.V,B,T,type = "put")
# B.t is the notional value of the debt discounted by the risk free rate at each timestep
# minus the BSM value of a put with spot prices V.t, rf = r, sig = sigma.V, strike = B.
lines(times491,Bt,col = "blue4") # Note how in the structural mode the value of the debt 
# converges to the total asset price as shareholder equity becomes worthless to cover 
# debtholder claims.


## Alternatively
set.seed(587)
V.t2 <- GBM(N = N, M = 1, x0 = V.0, t0 = 0, T = T, 
                theta = (mu.V - r)/sigma.V, sigma = sigma.V)
V.t2[365] # 1.47 > 0.85. No default, debtholders payed at T and
V.t2[365] - 0.85 # shareholder equity valued at 0.62
times587 <- seq(from = 0, to = 1, length = N+1)

plot(times587, V.t2, type = "l", ylim = range(0.6,max(V.t2)), xlab = "Timestep t",
     ylab = expression(V[t]))
abline(v = 1)
abline(h = .85)
# Default free debt value over t steps:
lines(times587, B*exp(-r*(T-times587)), col = "orange2")

B.t2 <-  B*exp(-r*(T-times587)) - Black_Scholes(t = times587, S = V.t2, r = r, 
                                                sigma = sigma.V, B, T, type = "put")
lines(times587,B.t2, col = "blue3")
## Note here value of debt converges to the default free ZCB price curve. 


# Referring to the text, Pr(V_T < B) = Pr(ln(V_T) < ln(B)):
pnorm((log(0.85/1)-(mu.V-0.5*sigma.V^2)*T)/sigma.V*sqrt(T))
## Under above assumptions and Merton coarse assumptions, 25.94% chance firm
# defaults at time T (that is iff net balance sheet price follows stochastic 
# diffusion model)


### Simulating Brownian Bridge ####

Bridge1 <- BB(N = 100, M = 1, x0 = 0, y = 0, t0 = 0, T = 1, Dt = NULL)

bridge.paths <- 500
bridge.matrix <- matrix(NA, nrow = 101, ncol = bridge.paths, byrow = FALSE)

for(i in 1:bridge.paths)
{
  bridge.matrix[,i] <- BB(N = 100, M = 1, x0 = 0, y = 0, t0 = 0, T = 1, Dt = NULL)
}

plot(bridge.matrix[,1], type = "l", ylim = c(-2,2), xlab = "Index", ylab = "Transition Level", 
     main = "Brownian Bridge: Equivalent Origination and Terminal Points (Simulated N = 500)")
for(i in 2:bridge.paths){
  lines(bridge.matrix[,i], type = "l", col = sample(1:20, 1))
}
