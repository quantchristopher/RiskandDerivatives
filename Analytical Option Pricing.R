##### Simple Analytical Options by Christopher Thomson


## The following code facilitates the pricing formulas for a variety of simple or "vanilla" financial options.
## Each contract is listed as a function of key inputs: 
## Underlying Asset Level  = S_t 
## Strike/Exercise Level   = K 
## Risk-Free Interest Rate = r
## Contract Maturity       = T
## Time to Maturity        = t


### Binary Option Solutions ####
## Call ption pays 1 unit of cash if in-the-money, nothing if out-of-the money
Binary.Call <- function(S_t,K,r,sigma,T,t) 
{
  d <- (log(S_t/K) + (r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  ## Compute
  exp(-r*(T-t))*pnorm(d, mean = 0, sd = 1)
}
# Example: A cash or nothing call with a strike of 58 on a stock thats currently trading
# out of the money at 55. Option expires in 233 units of time
Binary.Call(S_t = 55, K = 58, r = 0.03, sigma = 0.25, T = 1, t = 17*(1/250))

## Binary/Cash or Nothing Put:
Binary.Put <- function(S_t,K,r,sigma,T,t)
{
 d <- (log(S_t/K) + (r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
 # Compute:
 exp(-r*(T-t))*pnorm(-d, mean = 0, sd = 1)
}


### Analytical Solutions for options with a general cost of carry rate b: ####
### General structures for b:
## Non-dividend paying security                 b = r
## Continuous dividend paying yield q           b = r - q
## Foreign Currency with foreign risk free r_F  b = r - r_F
## Commodity Option with storage cost rate c_s  b = r + c_s
## Commodity Option with convenience yield y    b = r - y
## Option on a forward or futures contract      b = 0

## Vanilla Call Option with b inputs:
Euro.Call.b <- function(S_t, K, b, r, sigma, T, t)
{
  dplus <- (log(S_t/K)+(b+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  dminus <- (log(S_t/K)+(b-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  # Compute
  exp((b-r)*(T-t))*S_t*pnorm(dplus, mean = 0, sd = 1) - exp(-r*(T-t))*K*pnorm(dminus, mean = 0, sd = 1)
}

## Vanilla Put Option with b inputs:
Euro.Put.b <- function(S_t, K, b, r, sigma, T, t) 
{
  dplus <- (log(S_t/K)+(b+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  dminus <- (log(S_t/K)+(b-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  # Compute:
  exp(-r*(T-t))*K*pnorm(-dminus, mean = 0, sd = 1) - exp((b-r)*(T-t))*S_t*pnorm(-dplus, mean = 0, sd = 1)
}

## Example: Commodity Euro-Call with a store rate of 5%, 34 days into trading year:
Euro.Call.b(S_t = 92,K = 100, b = (.05 + .03), r = .03, sigma = 0.15, T = 1, t = 34/250)

## Example: Euro-Put on foreign currency with foreign risk free rate 2%, 45 days into trading year:
Euro.Put.b(S_t = 1, K = 1.55, b = (.03 - .02), r = .03, sigma = .11, T = 1, t = 45/250)


### Analytical Solution for stock that pays discrete dividend yield: ####
## Dividend yield is = q; number of dividend occurrences until option exercise/maturity = n
## b is arguably redudant here as it will just equal r for most assets

Euro.Call.discdiv <- function(S_t, K, b, r, q, n, sigma, T, t)
{
  dplus.disc <- (log((S_t*(1-q)^n)/K)+(b+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  dminus.disc <- (log((S_t*(1-q)^n)/K)+(b-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  
  exp((b-r)*(T-t))*(S_t*(1-q)^n)*pnorm(dplus.disc, mean = 0, sd = 1) - exp(-r*(T-t))*K*pnorm(dminus.disc, mean = 0, sd = 1)
}

## Note that for q = 0, n = 1 will give the normal Black-Scholes for a Euro-Call with
# the given underlying parameters:
Euro.Call.discdiv(S_t = 110, K = 100, b = .04, r = .04, q = 0, n = 1, sigma = 0.20, T = 0.5, t = 0)
# = 13.6955 ~ 13.70 per contract 

## Now compare to an equivalent option contract but the stock is known to pay a dividend yield of 2% in the next
# two and five months i.e. twice:
Euro.Call.discdiv(S_t = 110, K = 100, b = .04, r = .04, q = 0.02, n = 2, sigma = 0.20, T = 0.5, t = 0)
# = 10.33766 ~ 10.34 per contract

Euro.Put.divyield <- function(S_t, K, b, r, q, n, sigma, T, t)
{
  dplus.disc <- (log((S_t*(1-q)^n)/K)+(b+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  dminus.disc <- (log((S_t*(1-q)^n)/K)+(b-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
  # Compute:
  exp(-r*(T-t))*K*pnorm(dminus.disc, mean = 0, sd = 1) - exp((b-r)*(T-t))*(S_t*(1-q)^n)*pnorm(dplus.disc, mean = 0, sd = 1)
}


### Analytical Solution for Stock with fixed Cash Dividend with simple modified volatility ####
## The following inputs are required in addition to the common BSM inputs:
# D = Cash Dividend, N = Number of Dividends until option expiry, t_0 = Time of earliest next
# dividend, t_n = Time of final dividend before exercise/expiry. (Note that payments are assumed
# equally spaced over interval [t_0,t_n]. In practice the net discounted cash flow of future dividends
# may need to be computed for the actual payment schedule of the security). Its recommended to input 
# t_0 and t_n as fractions of the trading schedule's basis i.e. x/250. 
Euro.Call.cashdiv <- function(S_t, K, b, r, D, N, sigma, T, t, t_0, t_n)
{
  ti <- seq(from = t_0, to = t_n, by = (t_n-t_0)/N)
  disc.D <- rep(NA, N)
  for( i in 1:N) {
    disc.D[i] <- exp(-r*(t_n-ti[i]))*D
  }
  Dnet <- sum(disc.D) # Total NPV of future dividends 
  
  sigma.D <- sigma*(S_t/(S_t - Dnet)) # Simple Modified Volatility
  S_rt <- S_t - Dnet # Risky component of stock value
  
  dplus.cashdiv <- (log(S_rt/K)+(b+0.5*sigma.D^2)*(T-t))/(sigma.D*sqrt(T-t))
  dminus.cashdiv <- (log(S_rt/K)+(b-0.5*sigma.D^2)*(T-t))/(sigma.D*sqrt(T-t))
  
  exp((b-r)*(T-t))*(S_rt)*pnorm(dplus.cashdiv, mean = 0, sd = 1) - exp(-r*(T-t))*K*pnorm(dminus.cashdiv, mean = 0, sd = 1)
}

# Example: Say a stock is trading at 110 that pays a 1.60 dividend "roughly" every month. Consider Euro-Call with strike
# 100 that has 125 days (half-trading year) to maturity. Then there are 6 dividend payments due, say 6th one is due day 
# before the option can be exercised/expires. 
Euro.Call.cashdiv(110,100,0.02,0.02,1.60,6,.20,0.5,0,0,124/250)
# Cost of such a contract to an investor is 6.91.

Euro.Call.cashdiv(135,119,0.02,0.02,0.68,12,0.115,1,0,(21/250),(249/250))

### BSM Analytical Solution for Perpetual American Put ####
## This is effectively the only formulaic solution for an American Style derivative
# and it assumes several constants that may or may not exists empirically beyond a certain length of time
American.Put.perp <- function(S_t, K, r, sigma)
{
  ((sigma^2*K)/(2*r+sigma^2))*(K/(S_t*(1+(sigma^2/2*r))))^(2*r/sigma^2)
}

## Example: Consider a perpetual put on a stock that has absolutely tanked:
American.Put.perp(23, 85, .03, .14)
