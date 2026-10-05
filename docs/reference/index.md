# Package index

## Main fitting function

- [`HubbellGLM()`](https://alessandrozito.github.io/HubbellGLM/reference/HubbellGLM.md)
  : Hubbell regression This function mimics the GLM structure. We do so
  to avoid conflicts with the current glm implementation (if it gets
  updated)
- [`estimate_sigma()`](https://alessandrozito.github.io/HubbellGLM/reference/estimate_sigma.md)
  : Estimate the sigma parameter for a HubbellGLM model via BIC
  minimization

## Prediction

- [`predict_curve()`](https://alessandrozito.github.io/HubbellGLM/reference/predict_curve.md)
  : Predict the accumulation curve from a HubbellGLM model for new
  observations.

## Inference

- [`summary(`*`<HubbellGLM>`*`)`](https://alessandrozito.github.io/HubbellGLM/reference/summary.HubbellGLM.md)
  : Summary function for the HubbellGLM method
- [`vcov_shared()`](https://alessandrozito.github.io/HubbellGLM/reference/vcov_shared.md)
  : Variance-covariance matrix for a HubbellGLM fit, optionally adjusted
  for shared observations via a sandwich estimator.

## Family objects

- [`hubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/hubbell.md)
  : Hubbell family of generalized linear models
- [`quasihubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/quasihubbell.md)
  : Quasi-Hubbell family of generalized linear models

## Link and distribution functions

- [`make.link.hubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/make.link.hubbell.md)
  : Link function for the Hubbell generalized linear model
- [`inv_polyseries()`](https://alessandrozito.github.io/HubbellGLM/reference/inv_polyseries.md)
  : Function to invert the polyseries
- [`mean_dirichlet_process()`](https://alessandrozito.github.io/HubbellGLM/reference/mean_dirichlet_process.md)
  : Mean of the n. of species from Dirichlet process
- [`deriv_dirichlet_process()`](https://alessandrozito.github.io/HubbellGLM/reference/deriv_dirichlet_process.md)
  : Derivative of the mean of the Dirichlet process distribution
- [`inv_mean_dirichlet_process()`](https://alessandrozito.github.io/HubbellGLM/reference/inv_mean_dirichlet_process.md)
  : Invert the mean of the Dirichlet process distribution
- [`rDPspecies()`](https://alessandrozito.github.io/HubbellGLM/reference/rDPspecies.md)
  : Simulate the number of distinct species from a Dirichlet process

## Data

- [`BarroColorado`](https://alessandrozito.github.io/HubbellGLM/reference/BarroColorado.md)
  : Barro Colorado Island subplot data
