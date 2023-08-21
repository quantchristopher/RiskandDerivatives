##### Portfolio Optimization for given return and risk tolerances by Christopher Thomson

## The following code fulfills the task of a Markowitz portolio optimization problem for 
## n - risky assets. In other words, the goal is to efficiently allocate funds in a portfolio 
## to a number of assets with fully known return distributions and covariance/correlation relationships.
## By "efficiently" we mean to minimize the covariance/correlation risk of the portfolio while seeking an 
## a predetermined expected return on the entire portfolio. The only constraint on that must be maintained 
## is the "initial capital for investment" must be fully utilized and all available stocks must recieve allocation. 
## (however very very small and/or very very large allocations are allowed AND negative or shortselling is allowed). 
## The process is facilitated using Lagrangian Multipliers for matrix method solutions. 


### Function Setup #### 

## The function incorporates three variables: 1. The mean of returns vector - mu 2. The known covariance
## matrix - cv 3. The expected or "required" return on the aggregate portfolio - Er
markowitz <- function(mu,cv,Er) {
  n = length(mu)
  wuns = matrix(1,n,1)
  A = t(wuns)%*%solve(cv)%*%mu
  B = t(mu)%*%solve(cv)%*%mu
  C = t(wuns)%*%solve(cv)%*%wuns
  D = B*C - A^2
  lam1 = (C*Er-A)/D
  lam2 = (B-A*Er)/D
  wts = lam1[1]*(solve(cv)%*%mu) + lam2[1]*(solve(cv)%*%wuns)
  g = (B[1]*(solve(cv)%*%wuns)-A[1]*(solve(cv)%*%mu))/D[1]
  h = (C[1]*(solve(cv)%*%mu) - A[1]*(solve(cv)%*%wuns))/D[1]
  
  wts = g + h*Er
}

## An example of the above:
mu = matrix(c(0.02,0.10,0.20),3,1) 
n = length(mu)
cv = matrix(c(0.0001,0,0,0,0.04,0.02,0,0.02,0.16),n,n)
Er = 0.18
wts = markowitz(mu,cv,Er) # Store the output of the function
print(wts)
# Note the output of wts under these parameters requires that the first asset/mu[1] has
# to be sold short to minimized the variance.


### Tracing the efficient frontier of the portfolio for various return to risk tolerances by toggling E[rp]:
Er_vec <- matrix(seq(0.01,0.15,0.01),15,1)
sig_vec <- matrix(0,15,1)
j = 0
for (Er in Er_vec) {
  j = j+1
  wts = markowitz(mu,cv,Er)
  sig_vec[j] = sqrt(t(wts)%*%cv%*%wts)
}
# For given mu values of each component, the above code calculates the different weight
# allocations required to give an expected return on a portfolio and the associated standard
# deviation. 

## Plotting the efficient frontier:
plot(sig_vec, Er_vec, type = "l", main = "Efficient Frontier of Portfolios")
