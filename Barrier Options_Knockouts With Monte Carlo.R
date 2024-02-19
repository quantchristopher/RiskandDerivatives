##### Barrier Option Knockout Call via Monte Carlo Simulation by Christopher Thomson

### The following code provides a basic function for pricing a Knockout Barrier Call Option.
##  The contract's contingencies are relatively simple; it provides its holder European exercise rights 
##  to purchase a unit of stock for a given price on a given date. If on the given date the stock 
##  price is above/below a certain threshold then the option expires worthlessly. 

##  Note that this is a fairly crude pricing function as it relies on the averaging of an expected 
##  value of payoff discounted back to the present (or defined time instance). It also makes several specific 
##  assumptions about the behavior of the underlying asset, particularly that the asset price follows
##  a log-normal random walk with constant scaling factor of diffusion/randomness. This
##  factor is derived from the standard deviation of log-changes/returns over a previous period aka
##  the volatility of the asset. In reality volatility is a highly dynamic quantity that 
##  evolves rapidly. Likewise, the risk free interest rate of the option is also assumed to be static in this 
##  calculation when in reality the spot free risk interest rate is also a dynamic and constantly 
##  evolving quantity (albeit on that is considerably less random than volatility, arguably more deterministic).

### Function Parameters ####

# S_0 <- Initial Asset Level
# E   <- Exercise Price
# r   <- Risk Free Interest Rate
# sigma <- Asset Level Volatility
# Knockout <- Defined "Above or Below", breach of either triggers knockout of the option
# B   <- Barrier Level aka Asset Level Event Trigger
# T   <- Maturity
# t   <- Time to Maturity
# rebate <- Payment Received in Barrier Breach (For Knockout B == 0)
# M   <- Number of MC Simulations of Price Path

### Basic Function ####

Knockout.Option <- function(S_0, E, r, sigma, Knockout = c("Up","Down"), B, T, t, rebate, M)
{
  switch(Knockout,
         "Up" = {
           # Generate Final Stock Prices With M Simulations
           S_T = as.numeric(S_0*exp(rnorm(M, mean = (r-0.5*sigma^2)*(T-t), sd = sigma*sqrt(T-t))))
           
           # Initialize Payout Vector:
           Payoff.Vals = rep(0,M)
           
           # Perform MC Simulation
           for(i in 1:M){
             
             # Find where barrier is breached
             if(any(S_T[i] > B)) {
               Payoff.Vals[i] = rebate
               # If no barrier breach then compute 
             } else {
               Payoff.Vals[i] = max(S_T[i] - E, 0)
             }  # End If/Else Snippet   
             
           } # End For Loop 
           
           # Compute Average Value of Simulated Outcomes
           Option.Price <- exp(-r*(T-t))*(1/M)*sum(Payoff.Vals)
           
           return(Option.Price)
           
           
         }, # End "Call" for switch function  
         
         "Down" = {
           
           S_T = as.numeric(S_0*exp(rnorm(M, mean = (r-0.5*sigma^2)*(T-t), sd = sigma*sqrt(T-t))))
           
           
           Payoff.Vals = numeric(M)
           
           
           for(i in 1:M){
             
             
             if(any(S_T[i] < B)) {
               Payoff.Vals[i] = rebate
               
             } else {
               Payoff.Vals[i] = max(S_T[i] - E, 0)
             }  # End If/Else Snippet   
             
           } # End For Loop 
           
           # Compute Average Value of Simulated Outcomes
           Option.Price <- exp(-r*(T-t))*(1/M)*sum(Payoff.Vals)
           
           return(Option.Price)
           
         }, # End "Down" of switch function
         
         stop("Wrong 'Knockout'")
         
  ) # End switch function
  
} # End Function


## Note that when entering inputs, for "Up and Out" options, B > E needs to be maintained otherwise worthless asset
## For example:
Knockout.Option(S_0 = 110, E = 120, r = 0.06, sigma = 0.2, Knockout = "Up", B = 130, T = 1, t = (1/252), 
                rebate = 0, M = 100000)

## For "Down and Out", B < E needs to be maintained.
## For Example:
Knockout.Option(S_0 = 32, E = 29, r = 0.0666, sigma = 0.17, Knockout = "Down", B = 26, T = 1, t = (156/252), 
                rebate = 0, M = 100000)



## Plotting the price process visually:
X <- GBM(N = 252, M = 100, 
         x0 = 110, t0 = 0,  
         T = 1, theta = (0.06 - 0.5*0.2^2), 
         sigma = 0.2)

plot(X[,1], type = "l", ylim = c(30,200), xlab = "Time t", 
     ylab = "Asset Level S", main = "Discretized Evolution of S: 252 Observations per 100 Simulations")
for(i in 2:100){
  lines(X[,i], type = "l", col = sample(1:40,1))
}
abline(h = 130, lty = 2, lwd = 2, col = "red")
abline(h = 120, lwd = 2)
abline(v = 1.0)

legend(x = "topleft",
       legend = c("Strike E = 120","Barrier B = 130"),
       lty = c(1,2),
       lwd = c(2,2),
       col = c("black", "red"),
       cex = 0.66
)
