### Geomentric Brownian Motion Simulation by Christopher Thomson

## B-Motion can be calculated using the discrete form of the solution to the
## geometric brownian motion stochastic differential equation.
## S(t) is moved to S(t+h) where h is in years 

## S(t+h) = S(t)*e^[(mu-1/2sigma^2)h + sigma*L*sqrt(h)] where L~N(0,1) aka standard normal distribution RV

# Therefore let:
n <- 252 ## Number of trading days
S.0 <- 100 ## Spot price on day 1
mu <- 0.00 ## Spot Drift
sig <- 0.20 ## Spot Volatility
S <- matrix(0,1,n+1) ## Dummy matrix of 1 row, 252+1 columns all containing 0. NA would also be a sufficient entry
h <- (1/n) ## Daily time step on annualized basis

S[1] <- S.0 ## Entering intital price into dummy matrix

for (i in 2:(n+1)) {
  S[i]=S[i-1]*exp((mu-(0.5*sig^2)*h+sig*rnorm(1)*sqrt(h)))
}

plot(t(S)[,1], type = "l") # Plots GBM as a time series

# Load Package:
library(moments)

## Mapping descriptive statistics 
mean(s)
skewness(t(s))
kurtosis(t(s))

## 3 GBM random walks simulated in same lines of code:
O <- matrix(0,3,n+1) ## Creating feed-in matrix for multiple brownian motions
O[,1] <- S.0 ## Starting price is the same for all
for (j in seq(2,(n+1))) {
  O[,j]=O[,j-1]*exp((mu-0.5*sig^2)*h+matrix(rnorm(3),3,1)*sqrt(h))
}

## Plot 3 paths:
plot(t(O)[,1], ylim = c(0,200), type = "l")
lines(t(O)[,2], type = "l", col = "red")
lines(t(O)[,3], type = "l", col = "blue")

## Simulating 10000 paths with 252 "observances" of data:
B <- matrix(0,10000,n+1)
B[,1] <- S.0
for (j in seq(2,(n+1))) {
    B[,j] = B[,j-1]*exp((mu-0.5*sig^2)*h+matrix(rnorm(10000),10000,1)*sqrt(h))
}

## Terminal Values of every path:
finvec <- B[,253]
min(finvec) # Max
max(finvec) # Min

## Sort:
finvecsort <- sort(finvec)

# Distribution of terminal values:
Tabfin <- table(cut(finvecsort, breaks = 30))

## Plot of distribution:
plot(Tabfin)


