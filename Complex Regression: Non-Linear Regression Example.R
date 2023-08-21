##### Computing Forward-Rate Curve Using Non-Linear Regression:

## The following code utilizes treasury STRIP data to construct a forward-rate base yield curve. Non-linear regression is used.


### Empirically mapping STRIP prices against maturity as well as -ln(P(T_i))/delta T_i
TData <- read.table("strips_dec95.txt", header = TRUE)
STRIPData <- TData[order(TData$T),] ### Orders data by maturity from shortest to longest T
t <- seq(0,30, length = 100) ### Creates 100 entry vector increasing by .303 from 0 to 30 
empirical <- -diff(log(STRIPData$price))/diff(STRIPData$T)

plot(STRIPData$T, STRIPData$price, main = "Price vs. Maturity: Treasury STRIPS",
     xlab = "Maturity by Year", ylab = "STRIP ZCB Price", cex = 0.75) ### Price vs Maturity plot

plot(STRIPData$T[2:length(STRIPData$T)], empirical, ylim = c(0.025,0.075), lwd = 2,
    main = "EE Forward Rates", xlab = "Maturity in Years", 
    ylab = "Empirically Estimated Forward Rate", type = "b", cex = 0.75) ### Empirical Forward Rate vs Maturity plot

### Models to calculate D(T) using non linear regression function with discount prices

## r(t|theta) is an estimated quadratic function
fitQuad <- nls(STRIPData$price ~ 100*exp(-beta0*T - (beta1*T^2)/2 - (beta2*T^3)/3), 
               data = STRIPData, start = list(beta0 = 0.03, beta1 = 0, beta2 = 0) ) 

fitQuad2 <- nls(STRIPData$price ~ 100*exp(-beta0*STRIPData$T - (beta1*STRIPData$T^2)/2 - (beta2*STRIPData$T^3)/3), 
              data = STRIPdata, start = list(beta0 = 0.03, beta1 = 0, beta2 = 0) ) 
### Note that the pathway recognizes T from the data input. Also note the quadratic has
## already been integrated over. Also note that data input is required regardless of whether 
## STRIPData is added in the predictor input

Quadcoeffiecients <- summary(fitQuad)$coef[,1]
summary(fitQuad)
forwardrate_Quad <- Quadcoeffiecients[1] + (Quadcoeffiecients[2]*t) + (Quadcoeffiecients[3]*t^2)
### r(t|theta) is now a function given estimated betas from nls (listed as 100 values using t above
## as input)


## r(t|theta) is an estimated cubic function
fitcubic <- nls(STRIPData$price ~ 
                  100*exp(-beta0*T - (beta1*T^2)/2 - (beta2*T^3)/3 - (beta3*T^4)/4),
                data = STRIPData, start = list(beta0 = 0.03, beta1 = 0, beta2 = 0, beta3 = 0))
Cubiccoefficients <- summary(fitcubic)$coef[,1]
forwardrate_Cubic <- Cubiccoefficients[1] + (Cubiccoefficients[2]*t) + (Cubiccoefficients[3]*t^2) + (Cubiccoefficients[4]*t^3)

## r(t|theta) estimated as a quadratic spline function
fitSpline <- nls(STRIPData$price ~ 
                   100*exp(-beta0*T - (beta1*T^2)/2 - (beta2*T^3)/3 -
                             (T>15)*(beta3*(T-15)^3)/3 ), 
                 data = STRIPData, start = list(beta0 = 0.047, beta1 = 0.0024, beta2 = 0, beta3 = -0.00007) )
Splinecoefficients <- summary(fitSpline)$coef[,1]
forwardrate_Spline <- Splinecoefficients[1] + (Splinecoefficients[2]*t) + (Splinecoefficients[3]*t^2) + (t>15)*(Splinecoefficients[4]*(t-15)^2)


### Plot Comparing Quadratic vs Cubic vs Spline vs Empirical models of forward rate r(t)
plot(t, forwardrate_Quad, type = "l", ylim = c(0.025,0.075), lwd = 2, xlab = "Maturity", ylab = "Forward Rate")
lines(t, forwardrate_Cubic, lty = 2, lwd = 3, col = "red")
lines(t, forwardrate_Spline, lty = 3, lwd = 3, col = "blue")
points(STRIPData$T[2:length(STRIPData$T)], empirical, type = "b", cex = 0.5)

### Note that nls() may fail at higher order polynomial or exponential functions
## To estimate theta vector, it may be necessary to employ optim function using MLE

