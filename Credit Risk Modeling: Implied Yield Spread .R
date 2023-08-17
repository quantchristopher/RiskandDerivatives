##### Merton Implied Credit Spreads by Christopher Thomson

## The following code provides a simple approach to calculating the numerical rate spread between the 
## ccYTM of a debt issuing firm/company's defaultable debt over the ccYTM of a risk free debt security 
## with the "equivalent" term structure.

MertonSpread <- function(leverage, tau, sigma.V) 
{
  dt1 <- (-log(leverage) + 0.5*tau*sigma.V^2)/(sigma.V*sqrt(tau))

  dt2 <- dt1 - sigma.V*sqrt(tau)
  # Return the value:
  -log(pnorm(dt2) + (pnorm(-dt1))/leverage)/tau
}

## Assume fixed inputs of time to maturity and leverage; 
# volatility of assets ranges from low to high:
sig.V.range <- seq(from = 0.01, to = 0.5, length = 50)
cta <- MertonSpread(0.6, 2, sig.V.range)

plot(sig.V.range,cta,type = "l") # In short this demonstrates that as the volatility of
# returns on the total value of the assets {V_t} increases, the credit spread also
# increases exponentially. Creditors should demand higher YTM for riskier/more volatile
# issuers/obligors

## Now assume fixed leverage and fixed volatility; allow time to vary:
tau.V.range <- seq(from = 0.01, to = 5, length = 50)
ctb <- MertonSpread(leverage = 0.6, tau = tau.V.range, sigma.V = 0.25)
plot(tau.V.range, ctb, type = "l") # Credit spread/riskiness over risk free asset increases
# as time to maturity (of the theoretical ZCB) increases. 
### *** Note that this still does not capture intrinsic risks of short term credit spreads
# Reduced Form models are more useful when tau is small regardless of other inputs

## Finally assume fixed volatility and time to maturing; allow leverage to vary:
leverage.V.range = seq(from = 0.01, to = 0.99, length = 99)
# For low asset value volatility:
ctc.low <- MertonSpread(leverage = leverage.V.range, tau = 2, sigma.V = 0.15)
plot(leverage.V.range, ctc.low, type = "l") # Leverage has little effect on credit spread
# until it becomes significant proportion of financing for the firm even at low asset value 
# volatility
# For high asset value volatility:
ctc.high <- MertonSpread(leverage = leverage.V.range, tau = 2, sigma.V = 0.35)
plot(leverage.V.range, ctc.high, type = "l") # Smoother increase in credit spread as total 
# ratio of liabilities increase

ctc.high.long <- MertonSpread(leverage = leverage.V.range, tau = 5, sigma.V = 0.35)
plot(leverage.V.range, ctc.high.long, type = "l") # Earlier increase in credit spread,
# credit investor should therefore demand higher ccYTM for firms that issue a lot of debt
# for a long period of time, with high volatility in the continuous return on the firm's assets 
## *** Note in the book that short term credit spreads 
