##### Constructing Basic GARCH Models by Christopher Thomson

## The following code simulates simple ARMA + GARCH models 

### Simulating a basic AR(1) ARCH(1) Model ####
n <- 10200 # NOS
e <- rnorm(n) # 10200 RVs ~ N(0,1). When plotted as a time series this is the Gaussian White Noise Component
h <- e # Dummy Vector for h terms to be filled - once filled values will be the ARCH
y <- e # Dummy Vector
sigma.sqd <- e^2 # Dummy Vector to be filled as components for each h term of t 
omega <- 1 # Scaled unconditional autocovariance term of ARCH model aka gamma*h_infinity 
alpha <- 0.55 # Regressor coefficient for ARCH
phi <- 0.8 # Regressor coefficient for AR
mu <- 0.1 # Mean of the overt process

omega/(1-alpha) ; sqrt(omega/(1-alpha)) # Actual unconditional variance and standard deviation

for(t in 2:n)
{
  sigma.sqd[t] <- omega + alpha*h[t-1]^2 # Generates values of sigma^2 = omega + alpha*h_t-1^2
  h[t] <- e[t]*sqrt(sigma.sqd[t]) # Since h_t = e_t*sigma_t = e_t*sqrt(omega + alpha*h_t-1^2), calculates the next h
  y[t] <- mu + phi*(y[t-1] - mu) + h[t] # Coupled with AR(1) where phi is AR coefficient and h_t is innovation
}


plot(e[10001:n], type = "l", xlab = "t", ylab = expression(epsilon(t)), main = "Gaussian White Noise")
plot(sqrt(sigma.sqd[10001:n]), type = "l", xlab = "t", ylab = expression(sigma[t]), main = "Conditional Standard Deviation")
plot(h[10001:n], type = "l", xlab = "t", ylab = expression(h(t)), main = "ARCH(1)")
plot(y[10001:n], type = "l", xlab = "t", ylab = expression(y(t)), main = "AR(1) + ARCH(1)")

acf(h[10001:n], type = "correlation", plot = TRUE) 
acf(h[10001:n]^2, type = "correlation", plot = TRUE) # Here h(t)^2 has ACF akin to that of AR(1), where significant amount of autocorrelation
# occurs within the first lag
acf(y[10001:n], type = "correlation", plot = TRUE)



### Simulating basic ARCH(3) model ####
# **** Delete all stored variables from above before running this portion of script ****

n <- 10500 # NOS
e <- rnorm(n, mean = 0, sd = 1) # Gaussian White Noise with unit variance
h <- e # Numeric vector to be filled with ARCH(3) process values
sigma.sqd <- e # Numeric vector to be filled with sigma_t^2 value
alpha1 <- 0.097 # Volatility regression coefficient 1
alpha2 <- 0.031 # Volatility regression coefficient 2
alpha3 <- 0.002 # Volatility regression coefficient 3
omega  <- 0.988 


for(t in 4:n)
{
  sigma.sqd[t] <- omega + alpha1*(h[t-1]^2) + alpha2*(h[t-2]^2) + alpha3*(h[t-3]^2)
  h[t] <- e[t]*sqrt(sigma.sqd[t])
}

plot(e[10000:n], type = "l", main = "Gaussian White Noise") # Note e ~ N(0,1) for all t
plot(sqrt(sigma.sqd[10000:n]), type = "l", main = "Conditional Standard Deviation")
plot(h[10000:n], type = "l", main = "ARCH(3) Process")
lines(e[10000:n], type = "l", col = "red")


### Simulate AR(1) + GARCH(1,1) Process ####
# **** Delete all previously stored values ****
set.seed(722)

n <- 10500
e <- rnorm(n, mean = 0, sd = 1)
sigma.sqd <- e # Sigma^2 Process Values
h <- e # GARCH Process Values 
y <- e # AR + GARCH Process Values
omega <- 0.998 # "Future values of GARCH process" coefficient
alpha <- 0.091 # ARCH coefficient - Some contribution from previous values of {h_t} process
beta <- 0.886 # GARCH coefficient - Lots of contribution from previous values of sigma^2
phi <- 0.766 # Conditional Mean
mu <- 0.115 # Unconditional Mean

for(t in 2:n)
{
  sigma.sqd[t] <- omega  + alpha*h[t-1]^2 + beta*sigma.sqd[t-1]
  h[t] <- e[t]*sqrt(sigma.sqd[t])
  y[t] <- phi*(y[t-1] - mu) + h[t]
}

plot(e[10001:n], type = "l", main  = "Gaussian White Noise Process")
plot(sqrt(sigma.sqd[10001:n]), type = "l", main  = "Conditional Standard Deviation")
plot(h[10001:n], type = "l", main = "GARCH(1,1) Process")
plot(y[10001:n], type = "l", main = "AR(1) + GARCH(1,1) Process")
