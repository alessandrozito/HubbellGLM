# HubbellGLM

**HubbellGLM** implements Hubbell regression — a generalized linear
model (GLM) for $`\alpha`$-diversity. The model treats observed species
richness as a random variable governed by the Dirichlet process
(i.e. the Ewens sampling formula) and uses standard GLM machinery to
regress Hubbell’s fundamental biodiversity number $`\alpha`$ on
environmental covariates. A single fitted model yields coefficients with
standard errors and $`p`$-values, likelihood-based model selection
(AIC/BIC), and closed-form predictions of Shannon entropy, Simpson’s
index, and Hill numbers at any sample size.

## Installation

Install the development version from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("alessandrozito/HubbellGLM")
```

## Quick start

``` r

#-- Load the R packages
library(HubbellGLM)
library(vegan)
library(lmtest)

#-- Load the data
data(BCI)
data("BarroColorado")

# Estimate the growth-rate exponent sigma from the null model
formula <- "cbind(n, y) ~ EnvHet + Habitat + Stream"
sigma_hat <- estimate_sigma(formula, data = BarroColorado[1:40, ]) # -0.7852143
#sigma_hat <- -0.7852143

# Fit a covariate model
fit <- HubbellGLM(formula,
                  family = hubbell(sigma = sigma_hat),
                  data   = BarroColorado[1:40, ])
summary(fit)

# Predict richness at new sites
new_sites <- BarroColorado[41:50, ]
predict(fit, newdata = new_sites, type = "response")

# Species accumulation curve with confidence band
curve <- predict_curve(fit, xnew = new_sites, n = 600, npoints = 100)
plot(curve$n, curve$mean[, 1] + 1.96 * curve$se[, 1], lty = "dashed", type = "l")
lines(curve$n, curve$mean[, 1], col = "red")
lines(curve$n, curve$mean[, 1] - 1.96 * curve$se[, 1], lty = "dashed")

# Optional: two ways of adjusting standard errors.  
#=== Option 1 - Use the quasihubbell model
fitQuasi <- HubbellGLM(formula,
                  family = quasihubbell(sigma = sigma_hat),
                  data   = BarroColorado[1:40, ])
summary(fitQuasi)

#=== Option 2 - Adjust standard errors via shared species (Jaccard similarity)
brayDist <- vegdist(BCI[1:40, ])
Dshared <- 1 - as.matrix(2 * brayDist / (1 + brayDist))
vcov_fit <- vcov_shared(fit, similarity = Dshared)
coeftest(fit, vcov. = vcov_fit)
```

## Main functions

| Function | Description |
|:---|:---|
| [`HubbellGLM()`](https://alessandrozito.github.io/HubbellGLM/reference/HubbellGLM.md) | Fit a Hubbell regression (main fitting function) |
| [`hubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/hubbell.md) / [`quasihubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/quasihubbell.md) | GLM family objects |
| [`estimate_sigma()`](https://alessandrozito.github.io/HubbellGLM/reference/estimate_sigma.md) | Estimate the polynomial link exponent $`\sigma`$ |
| [`predict()`](https://rdrr.io/r/stats/predict.html) | Predict richness at new sites |
| [`predict_curve()`](https://alessandrozito.github.io/HubbellGLM/reference/predict_curve.md) | Species accumulation curve with standard errors |
| `vcov_species()` | Jaccard-adjusted sandwich variance–covariance matrix |

## Reference

Zito, A., Rigon, T., Roslin, T., Niittynen, P., Hebert, P. D. N.,
Zakharov, E. V., Ratnasingham, S., iBOL Consortium, Ovaskainen, O., and
Dunson, D. B. (2026). *Predicting global biodiversity via Hubbell
regression*. bioArxiv.
