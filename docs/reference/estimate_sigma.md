# Estimate the sigma parameter for a HubbellGLM model via BIC minimization

Estimate the sigma parameter for a HubbellGLM model via BIC minimization

## Usage

``` r
estimate_sigma(formula, data, startpoint = 0, verbose = TRUE)
```

## Arguments

- formula:

  A formula specifying the model (same syntax as
  [`HubbellGLM`](https://alessandrozito.github.io/HubbellGLM/reference/HubbellGLM.md)).

- data:

  A `data.frame` containing the variables in `formula`.

- verbose:

  Logical; if `TRUE` (default), prints the sigma value evaluated at each
  iteration.

## Value

A scalar: the value of `sigma` that minimizes the BIC of the fitted
model.
