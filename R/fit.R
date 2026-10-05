# Useful function that is present in the glm.fit...
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


# Function to fit the glm using the standard IRLS. Code is copied using the
# glm.fit function, adapting to include the parameter "size" and "sigma" in the code
glm.fit.hubbell <- function(x, y,
                            weights = rep.int(1, nobs),
                            start = NULL,
                            etastart = NULL,
                            mustart = NULL,
                            offset = rep.int(0, nobs),
                            family = hubbell(sigma = sigma),
                            control = list(),
                            intercept = TRUE,
                            singular.ok = TRUE) {
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) {
    rownames(y)
  } else {
    names(y)
  }
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
  }
  if (is.null(offset)) {
    offset <- rep.int(0, nobs)
  }
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)) {
    stop("'family' argument seems not to be a valid family object",
      call. = FALSE
    )
  }
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  valideta <- family$valideta %||% function(eta) TRUE
  validmu <- family$validmu %||% function(mu) TRUE
  # Check the argument sigma for the polynomial series
  if (!is.null(family$sigma)) {
    sigma <- family$sigma
  } else {
    sigma <- NULL
  }
  if (is.null(mustart)) {
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta)) {
      stop("invalid linear predictor values in empty model",
        call. = FALSE
      )
    }
    mu <- linkinv(eta, size, sigma)
    if (!validmu(mu)) {
      stop("invalid fitted means in empty model", call. = FALSE)
    }
    dev <- sum(dev.resids(y, mu, weights, size))
    w <- sqrt((weights * mu.eta(eta, size, sigma)^2) / variance(mu, size))
    residuals <- (y - mu) / mu.eta(eta, size, sigma)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  } else {
    coefold <- NULL
    eta <- etastart %||% {
      if (!is.null(start)) {
        if (length(start) != nvars) {
          stop(
            gettextf(
              "length of 'start' should equal %d and correspond to initial coefs for %s",
              nvars, paste(deparse(xnames), collapse = ", ")
            ),
            domain = NA
          )
        } else {
          coefold <- start
          offset + as.vector(if (NCOL(x) == 1L) {
            x * start
          } else {
            x %*% start
          })
        }
      } else {
        family$linkfun(mustart, size, sigma)
      }
    }
    mu <- linkinv(eta, size, sigma)
    if (!(validmu(mu) && valideta(eta))) {
      stop("cannot find valid starting values: please specify some",
        call. = FALSE
      )
    }
    devold <- sum(dev.resids(y, mu, weights, size))
    boundary <- conv <- FALSE
    for (iter in 1L:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu, size)[good]
      if (anyNA(varmu)) {
        stop("NAs in V(mu)")
      }
      if (any(varmu == 0)) {
        stop("0s in V(mu)")
      }
      mu.eta.val <- mu.eta(eta, size, sigma)
      if (any(is.na(mu.eta.val[good]))) {
        stop("NAs in d(mu)/d(eta)")
      }
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning(gettextf(
          "no observations informative at iteration %d",
          iter
        ), domain = NA)
        break
      }
      z <- (eta - offset)[good] + (y - mu)[good] / mu.eta.val[good]
      w <- sqrt((weights[good] * mu.eta.val[good]^2) / variance(mu, size)[good])
      fit <- .Call(stats:::C_Cdqrls, x[good, , drop = FALSE] *
        w, z * w, min(1e-07, control$epsilon / 1000), check = FALSE)

      if (any(!is.finite(fit$coefficients))) {
        conv <- FALSE
        warning(gettextf(
          "non-finite coefficients at iteration %d",
          iter
        ), domain = NA)
        break
      }
      if (nobs < fit$rank) {
        stop(sprintf(
          ngettext(
            nobs, "X matrix has rank %d, but only %d observation",
            "X matrix has rank %d, but only %d observations"
          ),
          fit$rank, nobs
        ), domain = NA)
      }
      if (!singular.ok && fit$rank < nvars) {
        stop("singular fit encountered")
      }
      start[fit$pivot] <- fit$coefficients
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset, size, sigma)
      dev <- sum(dev.resids(y, mu, weights, size))
      if (control$trace) {
        cat("Deviance = ", dev, " Iterations - ", iter,
          "\n",
          sep = ""
        )
      }
      boundary <- FALSE
      if (!is.finite(dev)) {
        if (is.null(coefold)) {
          stop("no valid set of coefficients has been found: please supply starting values",
            call. = FALSE
          )
        }
        warning("step size truncated due to divergence",
          call. = FALSE
        )
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit) {
            stop("inner loop 1; cannot correct step size",
              call. = FALSE
            )
          }
          ii <- ii + 1
          start <- (start + coefold) / 2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset, size, sigma)
          dev <- sum(dev.resids(y, mu, weights, size))
        }
        boundary <- TRUE
        if (control$trace) {
          cat("Step halved: new deviance = ", dev, "\n",
            sep = ""
          )
        }
      }
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold)) {
          stop("no valid set of coefficients has been found: please supply starting values",
            call. = FALSE
          )
        }
        warning("step size truncated: out of bounds",
          call. = FALSE
        )
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit) {
            stop("inner loop 2; cannot correct step size",
              call. = FALSE
            )
          }
          ii <- ii + 1
          start <- (start + coefold) / 2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset, size, sigma)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights, size))
        if (control$trace) {
          cat("Step halved: new deviance = ", dev, "\n",
            sep = ""
          )
        }
      }
      if (abs(dev - devold) / (0.1 + abs(dev)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      } else {
        devold <- dev
        coef <- coefold <- start
      }
    }
    if (!conv) {
      warning("glm.fit: algorithm did not converge", call. = FALSE)
    }
    if (boundary) {
      warning("glm.fit: algorithm stopped at boundary value",
        call. = FALSE
      )
    }
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps)) {
        warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
          call. = FALSE
        )
      }
    }
    if (family$family == "poisson") {
      if (any(mu < eps)) {
        warning("glm.fit: fitted rates numerically 0 occurred",
          call. = FALSE
        )
      }
    }
    if (fit$rank < nvars) {
      coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
    }
    xxnames <- xnames[fit$pivot]
    residuals <- (y - mu) / mu.eta(eta, size, sigma)
    fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
    } else {
      Rmat <- fit$qr[1L:nvars, 1L:nvars]
    }
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY) {
    names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int(
      "",
      sum(good) - fit$rank
    ))
  }

  wtdmu <- if (intercept) {
    # All models must have the same alpha here. This is estimated
    # through an appropriate function, which we specify now.
    #sum(weights * y) / sum(weights)
    beta_null <- estimateNullModel(y = y, n = size, sigma = sigma, family = family)
    linkinv(rep_len(beta_null, length(size)), size, sigma)
  } else {
    linkinv(offset, size, sigma)
  }
  nulldev <- sum(dev.resids(y, wtdmu, weights, size))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) {
    0
  } else {
    fit$rank
  }
  resdf <- n.ok - rank
  aic.model <- aic(y, n, mu, weights, dev, size) + 2 * rank
  list(
    coefficients = coef, residuals = residuals, fitted.values = mu,
    effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
    rank = rank, qr = if (!EMPTY) {
      structure(fit[c(
        "qr", "rank",
        "qraux", "pivot", "tol"
      )], class = "qr")
    }, family = family,
    linear.predictors = eta, deviance = dev, aic = aic.model,
    null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
    df.residual = resdf, df.null = nulldf, y = y, size = size, sigma = sigma,
    converged = conv,
    boundary = boundary, name_size = name_size, name_y = name_y
  )
}



# Fit the intercept-only model (M0) and return its coefficient beta0 = log(alpha).
#
# Every observation shares one alpha, but this is still not the usual glm
# shortcut sum(weights * y) / sum(weights): each observation has its own
# community size n_i, so the null means mu_i = g^-1(beta0; n_i) differ across
# observations even though beta0 does not.
#
# It is, however, nothing more than the ordinary IRLS of glm.fit.hubbell with
# X = 1, so the quantities are taken from the family object instead of being
# re-derived here. With a single column the weighted least squares step
# collapses to the scalar update beta0 <- sum(w * z) / sum(w), which is why no
# QR solve appears below.
estimateNullModel <- function(y, n, sigma, family, tol = 1e-10,
                              beta_start = NULL, maxiter = 100) {
  ny <- length(y)

  # eta_i solves mu(eta_i; n_i) = y_i one observation at a time. Since mu is
  # increasing in eta, beta0 = min(eta_i) makes every mu_i <= y_i and
  # beta0 = max(eta_i) makes every mu_i >= y_i, so the score
  #   U(beta0) = sum_i (y_i - mu_i) * (dmu_i/deta) / V(mu_i)
  # is non-negative at one end and non-positive at the other: the root is
  # bracketed by the per-observation solutions.
  eta_i <- family$linkfun(y, n, sigma)
  lo <- min(eta_i)
  hi <- max(eta_i)
  if (!(hi > lo)) return(mean(eta_i))

  beta <- if (is.null(beta_start)) mean(eta_i) else beta_start
  beta <- min(max(beta, lo), hi)

  for (t in seq_len(maxiter)) {
    eta <- rep_len(beta, ny)
    mu  <- family$linkinv(eta, n, sigma)
    dmu <- family$mu.eta(eta, n, sigma)      # d mu / d eta
    v   <- family$variance(mu, n)            # V(mu)

    w <- dmu^2 / v                           # IRLS weight
    z <- eta + (y - mu) / dmu                # working response
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) {
      stop("null model produced non-positive IRLS weights", call. = FALSE)
    }
    beta_new <- sum(w * z) / sw

    # beta_new - beta = U(beta) / sum(w), and sum(w) > 0, so the direction of
    # the IRLS step is the sign of the score. That narrows the bracket for free.
    if (is.finite(beta_new) && beta_new > beta) lo <- beta else hi <- beta
    # Unsafeguarded IRLS can be thrown a long way by a poor starting value
    # (mu near its boundary makes dmu tiny and z explode), so any step leaving
    # the bracket is replaced by a bisection.
    if (!is.finite(beta_new) || beta_new <= lo || beta_new >= hi) {
      beta_new <- 0.5 * (lo + hi)
    }

    converged <- abs(beta_new - beta) < tol * max(1, abs(beta_new))
    beta <- beta_new
    if (converged) break
  }
  if (t == maxiter) {
    warning("estimateNullModel: reached maxiter without converging",
            call. = FALSE)
  }
  beta
}

#' Variance-covariance matrix for a HubbellGLM fit, optionally adjusted for
#' shared observations via a sandwich estimator.
#'
#' @param fit An object of class \code{HubbellGLM}.
#' @param similarity An optional \eqn{n \times n} similarity matrix encoding
#'   dependence between the \eqn{n} observations used to fit \code{fit}.
#'   The diagonal must be 1 and all off-diagonal entries must be in
#'   \eqn{[0, 1)}. A sparse \code{Matrix} (\code{dgCMatrix} or the symmetric
#'   \code{dsCMatrix} returned by \code{get_shared_species}) is accepted and
#'   preferred at large \eqn{n}. If \code{NULL} (default), returns
#'   \code{vcov(fit)}.
#'
#' @return A \eqn{p \times p} variance-covariance matrix.
#'
#' @export
vcov_shared <- function(fit, similarity = NULL) {
  if (!inherits(fit, "HubbellGLM")) {
    stop("'fit' must be an object of class 'HubbellGLM'")
  }
  if (is.null(similarity)) {
    return(vcov(fit))
  }
  nobs <- nrow(model.matrix(fit))
  # A sparse Matrix is not a base matrix, so test for both.
  if (!(is.matrix(similarity) || methods::is(similarity, "Matrix"))) {
    stop("'similarity' must be a matrix or a Matrix")
  }
  if (!identical(dim(similarity), c(nobs, nobs))) {
    stop(sprintf(paste("'similarity' must be a %d x %d matrix matching the",
                       "number of observations in 'fit'"), nobs, nobs))
  }
  if (!isTRUE(all(Matrix::diag(similarity) == 1))) {
    stop("'similarity' must have 1 on the diagonal")
  }
  # Check the off-diagonal without ever forming it. The old
  # `similarity[row(similarity) != col(similarity)]` builds two n x n index
  # matrices and an n^2 - n vector, which at n = 15,700 is several GB spent
  # entirely on validation.
  if (methods::is(similarity, "sparseMatrix")) {
    # Zeros are structural and trivially in range, so only the stored values
    # need checking; the diagonal was already confirmed to be 1.
    st <- methods::as(similarity, "TsparseMatrix")
    off <- st@i != st@j
    bad_lo <- any(st@x[off] < 0)
    bad_hi <- any(st@x[off] >= 1)
  } else {
    rng <- range(similarity[upper.tri(similarity)],
                 similarity[lower.tri(similarity)])
    bad_lo <- rng[1] < 0
    bad_hi <- rng[2] >= 1
  }
  if (bad_lo || bad_hi) {
    stop("off-diagonal entries of 'similarity' must be in [0, 1)")
  }

  sigmaF      <- fit$sigma
  # The bread must be the UNSCALED inverse information (X'WX)^-1, not vcov().
  # For a quasi family vcov() is dispersion * (X'WX)^-1, and since the bread
  # appears on both sides it would contribute dispersion^2 while the meat below
  # - built from V(mu), not dispersion * V(mu) - contributes none. A sandwich
  # estimator is dispersion-free by construction; using vcov() here inflated
  # every standard error by exactly the dispersion (31.5 on the GMTP sample).
  FisherI_inv <- summary(fit)$cov.unscaled
  X           <- model.matrix(fit)
  y           <- fit$y
  mu          <- fit$fitted.values
  size        <- fit$size
  alpha       <- inv_mean_dirichlet_process(mu_target = mu, size = size)
  v           <- alpha * (digamma(alpha + size) - digamma(alpha)) +
                 alpha^2 * (trigamma(alpha + size) - trigamma(alpha))
  dlink       <- polyseries_meanvar_fast(size = size, alpha = c(exp(X %*% coef(fit))), sigma = sigmaF)$var
  # Row scaling, not a matrix product: diag(w) %*% X allocates an n x n dense
  # matrix (1.97 GB at n = 15,700) to do what recycling does for free.
  Scores      <- c((y - mu) / v * dlink) * X
  # `similarity` is only ever multiplied - it is never inverted, and never
  # needs to be, so a sparse operand stays sparse right through here.
  # Matrix::crossprod dispatches for a base matrix against a Matrix; base
  # crossprod() does not.
  meatD       <- as.matrix(Matrix::crossprod(Scores, similarity %*% Scores))
  FisherI_inv %*% meatD %*% FisherI_inv
}


