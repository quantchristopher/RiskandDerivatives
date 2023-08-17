##### Fitting Simple GARCH Models to Time Series Data by Christopher Thomson

## Fitting a simple GARCH model to BMW stock log-returns, built on simple 
## Gaussian White Noise Innovations.

# Needed dataset embedded in evir package --> install.packages("evir") 
library(evir) 
library(rugarch)
library(xts)

data(bmw, package = "evir")

BMW <- as.data.frame(bmw)
plot.ts(BMW, ylab = "Logarithmic Returns", main = "BMW Stock Returns")


ARMA.GARCH.normodel <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                                  variance.model = list(garchOrder = c(1,1))) 

# Above writes an AR(1) + GARCH(1,1) specified canvassing function to be applied to a dataset. 
# Here intended to apply to the BMW log-returns dataset.

BMW.ARGARCH.111 <- ugarchfit(data = BMW, spec = ARMA.GARCH.normodel) # Applies the spec model
show(BMW.ARGARCH.111)

# Having fit the GARCH model, can use the estimated parameters to simulate a stock price path

e <- rnorm(200, mean = 0, sd = 1)
h <- e
y <- e
alpha.1 <- 0.099405
beta.1 <- 0.863666
sigma.sqd <- e
omega <- 0.00009
phi <- 0.098133
mu <- 0.000453

for(t in 2:200)
{
  
}







