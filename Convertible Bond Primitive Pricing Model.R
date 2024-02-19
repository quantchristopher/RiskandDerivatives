##### Convertible Bond Primitive Pricing Model by Christopher Thomson

### Convertible bonds are hybrid financial instruments offering investors 
##  a complex debt instrument that can be converted to shares of equity if the 
##  owner so chooses to exercise. Accurate pricing of convertible bonds is complex
##  as they are effectively derivative contracts of two underlying state variables: 
##  the equity share price for which it may be converted, and the interest rate(s)
##  with which the future cash flows produced by the bond will be discounted by. 

## The R script below computes a primitive function for the price of a convertible bond
## in which the price is decomposed into the net present value of an equivalent non-convertible
## bond plus the net value of european call options on the underlying stock that would be
## needed to be exercised in order to buy the number of shares stated for the bond to convert to
## if conversion were to be exercised. Several simplifying assumptions are made and a 
## discussion on the implications of these assumptions is made at the bottom of this script


### Inputs ####
# S <- Current Stock Price
# r <- Continuous Risk Free Rate (can be assumed to be an averaged value over life of the bond) 
# sigma <- Stock Volatility (Annualized)
# T <- Maturity of the Underlying Bond
# t <- Assessed Time input
# n <- Conversion Ratio
# B <- Underlying Bond Redemption Value
# c <- Bond Coupon/Interest Payment
# y <- Subdivision of payments i.e. semiannual bond etc.
# q <- Credit Riskiness Spread Over Required over Risk Free Rate

CB.Primitive <- function(S, r, sigma, T, n, B, c, y, q){
  
  # Calculate cash flow of coupons (notional times rate over frequency of pmt)
  CF <- (B*c)/y
  
  # Calculate the value of the straight bond with as if it were not convertible
  Straight.Bond <- rep(0,T*y)
  Straight.Bond[T*y] <- (B + CF)*exp(-(r+q)*T)
  for(i in 1:((T*y)-1)){
    Straight.Bond[i] <- CF*exp(-(r+q)*(i*(1/y)))
  }
  # NPV
  SB <- sum(Straight.Bond)
  
  # Calculate "Strike" of Convertible Option:
  K <- (B/n)
  
  # Calculate value of European Call Option using BSM formula solution:
  dplus <- (log(S/K)+(r+0.5*sigma^2)*(T))/(sigma*sqrt(T))
  dminus <- (log(S/K)+(r-0.5*sigma^2)*(T))/(sigma*sqrt(T))
  CO <- S*pnorm(dplus, mean = 0, sd = 1) - exp(-r*T)*K*pnorm(dminus, mean = 0, sd = 1)
  
  # Take sum of the straight bond and n euro calls on to compute convertible bond price
  return(SB + n*CO)
  
  
}

## Test:
CB.Primitive(S = 47, r = 0.05, sigma = 0.22, T = 5, n = 20, B = 1000, c = 0.02, y = 2, q = 0)
# For reference - the above function should generate a value of 1125.53 


### While this model may stand as a quick solution to verifying the price of a simple 
##  convertible, its deficiencies become readily apparent when considering both the 
##  pricing function inputs and contingencies of the contract. First and foremost,
##  it is naive to assume that spot risk free interest rates will remain static over 
##  the course of 5 years; the actual price of the underlying bond is exposed to 
##  interest rate fluctuations that could severely impact the net present value of
##  its future cash flows. Likewise, the martingale discount rate of a call on the 
##  underlying stock is unlikely to remain a static single value. Additionally,
##  the volatility of the stock is assumed to be static when in reality this is a 
##  highly dynamic quantity. Likewise, the formulaic approach works best for long
##  maturities when input quantities can be averaged or concatenated over long periods.
##  There is also no input for embedded options on the bond which may drastically affect
##  the price and therefore effective yield.Possibly the most crucial aspect of valuation 
##  that is overlooked by this model is the underlying bond's credit quality. Presumably,
##  convertible bond issuers are likely to have some kind of credit risk with the prospect
##  of default. A "bandage" solution for this was to introduce a credit spread q over the
##  the risk free rate, however this fails to account for the correlative aspect the bond's
##  return performance will have with the associated stock price.



