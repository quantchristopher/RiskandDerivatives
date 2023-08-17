##### Simple Finite Differencing Function:

### The following code intends to replicate the explicit method of finite differencing
#   to solve for the intrinsic price of a derivative security with underlying S_t and
#   time to expiry t, i.e. V(S_t,t). It is assumed that V(S,t) is the solution to a 
#   PDE of the form:
#               (dV/dt) + a(S,t)*(d^2V/dS^2) + b(S,t)*(dV/dS) + c(S,t)*V(S,t) = 0. 

### Simple European Call/Put Options on underlying S ####
#   a(S,t) == 1/2*sigma^2*S^2 for everywhere t 
#   b(S,t) == r*S for everywhere t
#   c(S,t) == -r  for everywhere t

## Note that in order to properly iterate, S will be in terms of an asset step
#  step-size d.S multipled by an integer i from 1,...,I such that S.max = I*d.S
#  It is assumed that prices are log-normally distributed, therefore S.min = 0.
#  Likewise, t is in terms of d.t times integer k from 1,...,K such that expiry 
#  T - K*d.t = 0, and contract initialization T - 0*d.t = T

## To maintain stability of the numerical method, restrictions on the step sizes
#  of d.S and d.t must be set. For general BSM:
#  d.t <= 1/(sigma^2*I^2)
#  d.S <= 2a/|b|

### Euro-Call/Put as a function ####
## Following inputs will price a European Call or Put Option in the function Euro_FDPricer
# K     = Exercise Price
# r     = Risk Free Interest Rate (Deterministic Constant)
# sigma = Annualized Standard Deviation of Returns on the underlying S (Deterministic Constant)
# T     = Length of Contract
# NAS   = Number of Asset Steps - Will subdivide the output matrix's rows (ideally NAS > 500)


SxTau_LatticeSurface <- function(K, r, sigma, T, NAS, Type = c("Call","Put"))
{
  Type <- match.arg(Type)
  ## Inputs that are used for both calls and puts
  d.S <- (2*K)/NAS # An easy method for determining both S >> K and step size
  d.t <- 0.9/(sigma^2*NAS^2) # Bounded time slice 
  NTS <- round(T/d.t) + 1 # Number of time slices (integer required)
  d.t <- T/NTS # Exact time slice
  
  switch(Type,
         
         "Call" = { 
           # Dummy matrix to be filled
           Price.Grid <- matrix(data = NA, nrow = NAS, ncol = NTS, byrow = FALSE)
          
          # First for loop adds payoff/value of call option at T:
           for(i in 1:NAS)
           {
             Price.Grid[i,NTS] <- max((i*d.S - K), 0)
           } # End first initializing for loop
           
           for(k in NTS:1)
           { ## Two following lines establish boundary conditions:
             Price.Grid[1,k] <- 0 # For S = 0
             Price.Grid[NAS,k] <- (NAS*d.S - K)*exp(r*(k*d.t)) # For S = 2K >> K
             
             for(i in 2:(NAS-1)) 
             {
               a.Gamma <- ((0.5*sigma^2*(i*d.S)^2)*d.t)/(d.S^2)*(Price.Grid[(i+1),k] - 2*Price.Grid[i,k] + Price.Grid[(i-1),k])
               b.Delta <- ((r*(i*d.S)*d.t)/(2*d.S))*(Price.Grid[(i+1),k] - Price.Grid[(i-1),k])
               c.Rho   <- -r*d.t*Price.Grid[i,k]
               
               Price.Grid[i,(k-1)] <- Price.Grid[i,k] + a.Gamma + b.Delta + c.Rho
             } # End inner asset step for-loop
           } # End outer time slice for-loop
           
           Price.Grid  
           
         }, # End "Call" Type
         
         "Put" = {
           Price.Grid <- matrix(data = NA, nrow = NAS, ncol = NTS, byrow = FALSE)
           
           for(i in 1:NAS)
           {
             Price.Grid[i,NTS] <- max((K - i*d.S), 0)
           } # End first/initializing for-loop
           
           for(k in NTS:1)
           { ## Two following lines establish boundary conditions:
             Price.Grid[NAS,k] <- 0 # For S = 2K >> K
             Price.Grid[1,k] <- K*exp(r*(k*d.t)) # For S = 0
             
             for(i in 2:(NAS-1)) 
             {
               a.Gamma <- ((0.5*sigma^2*(i*d.S)^2)*d.t)/(d.S^2)*(Price.Grid[(i+1),k] - 2*Price.Grid[i,k] + Price.Grid[(i-1),k])
               b.Delta <- ((r*(i*d.S)*d.t)/(2*d.S))*(Price.Grid[(i+1),k] - Price.Grid[(i-1),k])
               c.Rho   <- -r*d.t*Price.Grid[i,k]
               
               Price.Grid[i,(k-1)] <- Price.Grid[i,k] + a.Gamma + b.Delta + c.Rho
             } # End inner asset step for-loop
           } # End outer time slice for-loop
           
           Price.Grid # Print Result
         },
         stop("Wrong 'Type' Input")
  ) # End Switch
} # End Function

# *** For simplicity in displaying results its ideal to set the function as a variable i.e. X <- Euro_FDPricer(...)









LuckyStrike_25 <- SxTau_LatticeSurface(25,0.03,0.11,1,500,"Call")
