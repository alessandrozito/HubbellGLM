# Hubbell family of generalized linear models

Hubbell family of generalized linear models

## Usage

``` r
hubbell(link = "polyseries", sigma = 0)
```

## Arguments

- link:

  The link function to use. Available options are `dp`, and `polyseries`

- sigma:

  Hyperparameter determining the polynomial growth. Must be `sigma < 1`
  for link `link = 'polyseries'`. It is set to 0 when `link = 'dp'`
  (logarithmic growth).
