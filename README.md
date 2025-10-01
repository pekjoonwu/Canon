# Canon (CAsual relationship identificatioN using ONe sample instrumental variable model)
## Installation
To install the latest version of the Canon package from GitHub, run the following code in R:
```
install.packages('devtools')
library(devtools)
devtools::install_github('pekjoonwu/Canon')
```
This command should automatically install all required packages if they are not installed already.

## Quick Start
```
library(Canon)
?run_Canon
```

## Usage
```
###
## Simulate gRNA matrix, exposure gene expression and outcome gene expression
###
library(MASS)
set.seed(2025)
n <- 10000
p <- 100

## Simulate the genetic components
simulate_matrix <- function(n, p, max_ones = 3) {
  set.seed(2025)
  mat <- matrix(0, n, p)
  mat <- t(apply(mat, 1, function(x) {
    k <- sample(0:max_ones, 1)
    if (k > 0) x[sample(seq_along(x), k)] <- 1
    x
  }))
  mat
}

grna_mat <- simulate_matrix(n, p)
beta <- rep(0, p)
true_idx <- sample(1:p, round(p * 0.3), replace = F)
beta[true_idx] <- rnorm(length(true_idx), 0, sqrt(0.025/length(true_idx)))

## calcualte the alpha value
alpha <- 0.3

## Simulate the random errors and data
sigma <- matrix(c(1 - 0.025, 0.15 * (sqrt(1 - 0.025) * sqrt(1 - 1e-3)), 
          0.15 * (sqrt(1 - 0.025) * sqrt(1 - 1e-3)), 1 - 1e-3), 2, 2)
dat <- do.call(rbind, lapply(1:n, function(x) {
  mu <- matrix(c(grna_mat[x, ] %*% beta, grna_mat[x, ] %*% beta * alpha), 2, 1)
  bvn1 <- mvrnorm(1, mu = mu, Sigma = sigma)
  return(bvn1)
}))

x <- dat[, 1]
y <- dat[, 2]

###
## Run Canon
###
results <- run_Canon(x_expression = as.matrix(x), 
                        y_expression = as.matrix(y), gRNA = grna_mat,
                        Gibbsnumber = 2000)
```
