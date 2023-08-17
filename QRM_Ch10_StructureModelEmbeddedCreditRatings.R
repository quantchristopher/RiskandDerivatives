##### QRM Chapter 10 - Embedding Credit Rating Migration into Structural Model

#### Packages and Setup ####
library(expm)

## Compute Generator Matrix:
Lambda <- matrix(c(0,       0.082,  0.0063, 0,      0.0003, 0,      0,      0,      0,
                   0.0091,  0,      0.0843, 0.0049, 0.0006, 0.0002, 0.0001, 0,      0.0002,
                   0.0006,  0.0248, 0,      0.0547, 0.0057, 0.0011, 0.0003, 0,      0.0006,
                   0.00039, 0.0017, 0.0411, 0,      0.0405, 0.0755, 0.0163, 0.0002, 0.0017,
                   0.0001,  0.0005, 0.0035, 0.0552, 0,      0.0722, 0.0058, 0.0007, 0.0106,
                   0.0001,  0.0003, 0.0011, 0.0032, 0.0458, 0,      0.0581, 0.0059, 0.0385,
                   0.0001,  0.0002, 0.0002, 0.0012, 0.0038, 0.0870, 0,      0.0372, 0.1334,
                   0,       0,      0,      0,      0.004,  0.0203, 0.0938, 0,      0.3793,
                   0,       0,      0,      0,      0,      0,      0,      0,      0),
                 ncol = 9, byrow = TRUE)

ratings <- c("Aaa","Aa","A","Baa","Ba","B","Caa","C","D")
dimnames(Lambda) <- list(ratings,ratings)

Lambda

Lambda.d <- rowSums(Lambda) # Sums values by row in "Aaa" "Aa" ... etc. 
diag(Lambda) <- Lambda.d # Above Lambda matrix has zero diagonals 
#  but insert rowsums as diagonal values into Lambda
rowSums(Lambda) 

P <- expm(Lambda) # Exponential function raised to power of matrix
rowSums(P)

P2 <- P[nrow(Lambda):1,ncol(Lambda):1] # computes reverse orders of rows/columns of P
P2 <- P2[-1,] # Eliminates first row

## Pick bond rated C:
(C.probs <- P2[1,]) # Transition probabilites for a firm rated at C initially
(C.cumulative <- cumsum(Cprobs))
C.thresholds <- vector(mode = "numeric", length = 9L) 

for(i in  1:length(C.cumulative))
{
  qnorm_prob <- qnorm(C.cumulative[i])
  C.thresholds[i] <- qnorm_prob
}


0.941*(-2.5) + 2.380952

.024/1.05

.029875/1.05
