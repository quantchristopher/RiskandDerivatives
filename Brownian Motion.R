### Brownian Motion Practice ###

## B-Motion can be calculated using the discrete form of the solution to the
## geometric brownian motion stochastic differential equation.
## S(t) is moved to S(t+h) where h is in years 

## S(t+h) = S(t)*e^[(mu-1/2sigma^2)h + sigma*L*sqrt(h)] where L~N(0,1) aka standard normal distribution RV

## Let:

n <- 252 ## Number of trading days
so <- 100 ## Spot price on day 1
mu <- 0.00 ## Spot Drift
sig <- 0.20 ## Spot Volatility
s <- matrix(0,1,n+1) ## Dummy matrix of 1 row, 252+1 columns all containing 0
h <- (1/n) ## Daily time step on annualized basis

s[1] <- so ## Entering intital price into dummy matrix

for (i in 2:(n+1)) {
  s[i]=s[i-1]*exp((mu-(0.5*sig^2)*h+sig*rnorm(1)*sqrt(h)))
}

plot(t(s)[,1], type = "l")

library(moments)

## Mapping descriptive statistics 
mean(s)
skewness(t(s))
kurtosis(t(s))

o <- matrix(0,3,n+1) ## Creating feed-in matrix for multiple brownian motions
o[,1] <- so ## Starting price is the same for all
for (j in seq(2,(n+1))) {
  o[,j]=o[,j-1]*exp((mu-0.5*sig^2)*h+matrix(rnorm(3),3,1)*sqrt(h))
}

plot(t(o)[,1], ylim = c(0,200), type = "l")
lines(t(o)[,2], type = "l", col = "red")
lines(t(o)[,3], type = "l", col = "blue")

b <- matrix(0,10000,n+1)
b[,1] <- so
for (j in seq(2,(n+1))) {
    b[,j] = b[,j-1]*exp((mu-0.5*sig^2)*h+matrix(rnorm(10000),10000,1)*sqrt(h))
}

finvec <- b[,253]
min(finvec)
max(finvec)

finvecsort <- sort(finvec)

Tabfin <- table(cut(finvecsort, breaks = 30))

plot(Tabfin)






## Converting independent random variables to bivariate random variable use following transform
## If random variables (e1,e2) ~ N(0,1) are independent, convert them to (x1,x2) with cor rho
## such that x1 = e1, x2 = rho*e1 + sqrt(1-rho^2)*e2

## Generate 10000 pairs of correlated random variates using following:

e <- matrix(rnorm(20000),10000,2)
cor(e)

cor(e[,1],e[,2])
