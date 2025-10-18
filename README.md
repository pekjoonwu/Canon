# Canon (CAsual relationship identificatioN using ONe sample instrumental variable model)
![scheme](https://github.com/pekjoonwu/Canon/blob/master/overview.png)

<p align="justify">

A critical analytical task in sc-CRISPR screening is identifying downstream genes influenced by perturbed target genes. Existing methods for this task primarily rely on traditional association- based analyses, which not only fall short in establishing causal relationships but also suffer from high false positive rates and limited statistical power. To overcome these limitations, we introduce a novel causal inference based framework that leverages the perturbation status of gRNAs in single cells as instrumental variables (IVs) to infer causal gene relationship via IV analysis. Building upon this framework, we further present Canon, a one-sample IV analysis method specifically tailored to systematically identify genes that are potentially causally influenced by perturbed target genes across diverse sc-CRISPR platforms. Canon ensures robust type I error control while maintaining high statistical power. We evaluated its performance through comprehensive simulations and real data applications. The gene-gene relationships identified by Canon provide valuable insights into the causal gene regulatory network, uncovering novel therapeutic targets for cancer treatment and demonstrating the transformative potential of sc-CRISPR screening to resolve causal networks at an unprecedented scale.

</p>

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

## Detailed Tutorial
See [Tutorial](https://pekjoonwu.github.io/Canon-analysis/) for detailed documentation and examples.


## Simple Usage Example
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

grna_mat <- scale(simulate_matrix(n, p))
beta <- rep(0, p)
true_idx <- sample(1:p, round(p * 0.3), replace = F)
beta[true_idx] <- rnorm(length(true_idx), 0, sqrt(0.025/length(true_idx)))
beta_off <- rep(0, p)
beta_off[true_idx] <- rnorm(length(true_idx), 0, sqrt(0.005/length(true_idx)))


## calcualte the alpha value
alpha <- sqrt(1e-3/0.025)

## Simulate the random errors and data
sigma <- matrix(c(1 - 0.025, 0.15 * (sqrt(1 - 0.025) * sqrt(1 - 1e-3 - 0.005)), 
          0.15 * (sqrt(1 - 0.025) * sqrt(1 - 1e-3 - 0.005)), 1 - 1e-3), 2, 2)
dat <- do.call(rbind, lapply(1:n, function(x) {
  mu <- matrix(c(grna_mat[x, ] %*% beta, grna_mat[x, ] %*% beta * alpha + grna_mat[x, ] %*% beta_off), 2, 1)
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
                        Gibbsnumber = 1000)
```
