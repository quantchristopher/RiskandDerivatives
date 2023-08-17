##### Portfolio Optimization for given return and risk tolerances

## The following code fulfills the task of a Markowitz portolio optimization for 
## n - risky assets using Lagrangian Multipliers and matrix method solutions. 


### Function Setup #### 
## Markowitz portfolio optimization of w_vec, given constraint that all of
## w_vec must be invested, portfolio covariance matrix must be minimized

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

mu = matrix(c(0.02,0.10,0.20),3,1) 
n = length(mu)
cv = matrix(c(0.0001,0,0,0,0.04,0.02,0,0.02,0.16),n,n)
Er = 0.18

wts = markowitz(mu,cv,Er)
print(wts)
### Note the output of wts under these parameters requires that the first asset/mu[1] has
## to be sold short to minimized the variance

### Quadprog Package can also be used to conduct these optimization problems with the added
## benefit of including additional constraints such as inequalities such as no-shortselling,
## minimum variance must still be above a certain threshold etc.


### Tracing the efficient frontier by toggling E[rp]
Er_vec <- matrix(seq(0.01,0.15,0.01),15,1)
sig_vec <- matrix(0,15,1)
j = 0
for (Er in Er_vec) {
  j = j+1
  wts = markowitz(mu,cv,Er)
  sig_vec[j] = sqrt(t(wts)%*%cv%*%wts)
}
### For given mu values of each component, the above code calculates the different weight
## allocations required to give an expected return on a portfolio and the associated standard
## deviation. 
 
plot(sig_vec, Er_vec, type = "l", main = "Efficient Frontier of Portfolios")
