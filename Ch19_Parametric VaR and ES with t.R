##### Chapter 19 - VaR and ES using parametric fitting 

## Recall: let F(r|thetavec) be a parametric family of distributions used to model the return
## distribution given parameters thetavec. If thetavec is estimated via MLE from historic returns
## then F^-1(alpha|thetavec) is an estimate of the alpha-quantile of the return distribution
## and VaR_parest(alpha) = -S * F^-1(alpha|thetavec) where S is the size of the current position 
## before the holding period of the return distribution return variabe will take place. 

## Also recall: if f(r|thetavec) is the density of the distribution -> F(r|thetavec) 
## then the ES_parest(alpha) = -(S/alpha) * integral[-infinity,F^-1(alpha|thetavec)](x*f(x|thetavec)dx)

### Example 19.3 - Assume that returns are iid t-distributed such that
## VaR(alpha) = -S * [mu_hat + q_alpha/t(nu_hat)*lambda_hat]
## which mu = mean, nu = scale, lambda = tail index of a sample of returns 
## and q is the alphath quantile of the t-distribution with index nu so that the 
## entire [""] quantity is the alpha-quantile of the loss distribution

data(SP500, package = "Ecdat")
### 2783 observations of the SP500 daily return
n <- 2783
SPreturn <- SP500$r500[(n-999):n] ### Shortens data to the last 1000 returns
year = 1981 + (1:n)*((1991.25-1981)/n) ### Creates sequence of 2783 units increasing by 0.00368307582
## from 1981 onward

alpha = 0.05
library(MASS)

fit_t <- fitdistr(SPreturn, "t")
param <- as.numeric(fit_t$estimate) ### Stores mu, sigma, and nu as numerical valued vector
mu <- param[1]
nu <- param[3] ### Degrees of Freedom  
lambda <- param[2]
sd <- param[2]*sqrt((nu)/(nu-2))
### Note that in computing theta_hat here; only mu, nu, and lambda are given

qalpha <- qt(alpha, df = nu) ### mu + lambda*qalpha is effectively F^-1(alpha|theta) 
## under parameters estimated from above

### Assuming an inital position of $20,000 in an SP500 Index

VaR_par_t <- -20000*(mu + lambda*qalpha) ### Indicates $324.17 is at the 5% CI risk
### In other words, theres a 5% chance the investment in SP500 will lose $324.17 or more
## over the holding period 

### For Estimated Shortfall
es1 <- dt(qalpha, df = nu) / alpha
es2 <- (nu+qalpha^2) / (nu - 1)
es3 <- -1*mu + lambda *es1*es2 

ES_par <- 20000*es3 ### There is a 5% chance the position will lose $324.17 or more
## The expected loss of within this 5% confidence range is $543.81 


###### Computing VaR and ES of univariate distribution of returns using empirical data
### with conditional mean and conditional variance adjustments i.e. ARMA + GARCH

library(rugarch)
garch.t <- ugarchspec(mean.model = list(armaOrder = c(1,0)), ## AR(1)
                      variance.model = list(garchOrder = c(1,1)), ## GARCH(1,1)
                      distribution.model = "std") ### Specifies ARMA + GARCH Model, empty values


sp.garch.t <- ugarchfit(data = SPreturn, spec = garch.t) ### Fits above model to given data, 
## estimating mu, phi, omega, alpha, beta, and shape nu ---> most interested in nu
show(sp.garch.t)

pred <- ugarchforecast(sp.garch.t, data = SPreturn, n.ahead = 5)
pred

nu2 <- as.numeric(coef(sp.garch.t)[6]) ### Collects shape parameter nu from GARCH model
q <- qstd(alpha, mean = 0.0007744, sd = 0.009486, nu = nu2) 
### Above computes the .05 quantile of the conditional distribution of the next return, estimated from 
## the predictive garch model

### Assuming a $20k initial investment in the SP500:
VaR_arma_garch <- -20000*q
### At this 5% confidence interval, theres a 5% chance the investment in the SP500 will
## lose $275.40 or more over the 1 period hold i.e. from t to t+1

### For Estimated Shortfall
lambda2 <- 0.009486/sqrt((nu2)/(nu2-2))
qalpha2 <- qt(alpha, df = nu2)
es1_2 <- dt(qalpha, df = nu2)/(alpha)
es2_2 <- (nu2 + qalpha2^2)/(nu2-1)
es3_2 <- -0.0007744 +lambda2*es1_2*es2_2

ES_par2 <- 20000*es3_2

dim(sigma(sp.garch.t)) ### 1000 estimates of sigma based on the arma+garch t-model
sig_avg <- mean(sigma(sp.garch.t))

x <- seq(1,1000, 1)
plot(x, sigma(sp.garch.t), type = "l")





