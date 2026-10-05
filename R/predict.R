#' Function to predict the response from a HubbellGLM model.
#' @exportS3Method
#' @keywords internal
predict.HubbellGLM <- function(
    object, newdata = NULL, type = c(
      "link", "response",
      "terms"
    ), se.fit = FALSE, dispersion = NULL, terms = NULL,
    na.action = na.pass, ...) {
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL
  if (!se.fit) {
    if (missing(newdata)) {
      pred <- switch(type,
        link = object$linear.predictors,
        response = object$fitted.values,
        terms = predict.lm(object,
          se.fit = se.fit, scale = 1, type = "terms",
          terms = terms
        )
      )
      if (!is.null(na.act)) {
        pred <- napredict(na.act, pred)
      }
    } else {
      pred <- predict.lm(object, newdata, se.fit,
        scale = 1,
        type = if (type == "link") {
          "response"
        } else {
          type
        }, terms = terms, na.action = na.action
      )
      switch(type,
        response = {
          if (object$family$family != "hubbell") {
            pred <- family(object)$linkinv(pred)
          } else {
            sigma <- family(object)$sigma
            pred <- family(object)$linkinv(pred, newdata[[object$name_size]], sigma)
          }
        },
        link = ,
        terms =
        )
    }
  } else {
    if (inherits(object, "survreg")) {
      dispersion <- 1
    }
    if (is.null(dispersion) || dispersion == 0) {
      dispersion <- summary(object, dispersion = dispersion)$dispersion
    }
    residual.scale <- as.vector(sqrt(dispersion))
    pred <- predict.lm(object, newdata, se.fit,
      scale = residual.scale,
      type = if (type == "link") {
        "response"
      } else {
        type
      }, terms = terms, na.action = na.action
    )
    fit <- pred$fit
    se.fit <- pred$se.fit
    switch(type,
      response = {
        if (object$family$family != "hubbell") {
          se.fit <- se.fit * abs(family(object)$mu.eta(fit))
          fit <- family(object)$linkinv(fit)
        } else {
          sigma <- family(object)$sigma
          if (missing(newdata)) {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit, object$size, sigma))
            fit <- family(object)$linkinv(fit, object$size, sigma)
          } else {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit, newdata[[object$name_size]], sigma))
            fit <- family(object)$linkinv(fit, newdata[[object$name_size]], sigma)
          }
        }
      },
      link = ,
      terms =
      )
    if (missing(newdata) && !is.null(na.act)) {
      fit <- napredict(na.act, fit)
      se.fit <- napredict(na.act, se.fit)
    }
    pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
  }
  pred
}

#' Predict the accumulation curve from a HubbellGLM model for new observations.
#'
#' @param fit An object of class \code{HubbellGLM}.
#' @param xnew A \code{data.frame} with the same predictor columns used to fit
#'   \code{fit}. Each row produces one set of accumulation curve predictions.
#' @param n Total number of individuals to extrapolate to. Must be a single
#'   positive integer.
#' @param npoints Number of equally-spaced grid points between 1 and \code{n}.
#'   Must satisfy \code{npoints <= n}. Default is 100.
#' @param .vcov Optional variance-covariance matrix for the regression
#'   coefficients. Must be a square matrix of dimension equal to the number of
#'   coefficients in \code{fit}. Defaults to \code{vcov(fit)}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{n}}{Integer vector of length \code{npoints} with the grid
#'     points (number of individuals).}
#'   \item{\code{mean}}{A \code{data.frame} with \code{npoints} rows and one
#'     column per row of \code{xnew}, containing the predicted species
#'     richness at each grid point.}
#'   \item{\code{se}}{A \code{data.frame} of the same shape as \code{mean},
#'     containing standard errors computed via the delta method.}
#' }
#'
#' @export
predict_curve <- function(fit, xnew, n = 3000, npoints = 100, .vcov = NULL) {
  if (!inherits(fit, "HubbellGLM")) {
    stop("'fit' must be an object of class 'HubbellGLM'")
  }
  if (!is.data.frame(xnew)) {
    stop("'xnew' must be a data.frame")
  }
  if (!is.numeric(n) || length(n) != 1L || n < 1L) {
    stop("'n' must be a single positive integer")
  }
  if (!is.numeric(npoints) || length(npoints) != 1L || npoints < 1L) {
    stop("'npoints' must be a single positive integer")
  }
  if (npoints > n) {
    stop("Must have npoints <= n")
  }
  if (!is.null(.vcov)) {
    p <- length(coef(fit))
    if (!is.matrix(.vcov) || !identical(dim(.vcov), c(p, p))) {
      stop(sprintf("'.vcov' must be a %d x %d matrix matching the number of coefficients in 'fit'", p, p))
    }
  }

  # Grid of points of the accumulation curve
  points_grid <- round(seq(from = 1, to = n, length.out = npoints)) - 1 # Adjust for zero due to Cpp

  # Full model matrix: nrow(xnew) x p
  Xmat <- model.matrix(terms(fit), data = xnew, xlev = fit$xlevels)
  # Variance-covariance matrix
  vcov_pred <- if (is.null(.vcov)) vcov(fit) else .vcov

  betas    <- coef(fit)
  nobs     <- nrow(xnew)
  mean_mat <- matrix(NA_real_, nrow = npoints, ncol = nobs)
  se_mat   <- matrix(NA_real_, nrow = npoints, ncol = nobs)

  for (i in seq_len(nobs)) {
    Xpred_i   <- Xmat[i, ]
    alpha_i   <- c(exp(betas %*% Xpred_i))
    se_fit_i  <- sqrt(c(Xpred_i %*% vcov_pred %*% Xpred_i))
    g_prime_i <- polyseries_var_grid(points_grid, size = n, alpha = alpha_i, sigma = fit$sigma)
    mean_mat[, i] <- polyseries_mean_grid(points_grid, size = n, alpha = alpha_i, sigma = fit$sigma)
    se_mat[, i]   <- g_prime_i * se_fit_i
  }

  colnames(mean_mat) <- rownames(Xmat)
  colnames(se_mat)   <- rownames(Xmat)

  return(list(
    n    = points_grid + 1L,
    mean = as.data.frame(mean_mat),
    se   = as.data.frame(se_mat)
  ))
}



# predict_curve <- function(fit, xnew, n = 3000, npoints = 100, .vcov = NULL) {
#   if(npoints > n) {
#     stop("Must have npoints <= n")
#   }
#   # Grid of points of the accumulation curve
#   points_grid <- round(seq(from = 1, to = n, length.out = npoints)) - 1 # Adjust for zero due to Cpp
#
#   # Find predictor in fit format
#   Xpred <- c(model.matrix(terms(fit), data = xnew, xlev  = fit$xlevels))
#   # Variance-covariance matrix
#   if(is.null(.vcov)) {
#     vcov_pred <- vcov(fit)
#   } else {
#     vcov_pred <- .vcov
#   }
#   # Obtain standard error for linear predictor
#   betas <- coef(fit)
#   alpha <- c(exp(betas %*% Xpred))
#   se_fit <- sqrt(c(Xpred %*% vcov_pred %*% Xpred))
#   # Adjust factor (delta method)
#   g_prime <- polyseries_var_grid(points_grid, size = n, alpha = alpha, sigma = fit$sigma)
#   # Caculate the prediction curve at the grid points
#   pred <- polyseries_mean_grid(points_grid, size = n, alpha = alpha, sigma = fit$sigma)
#   # Return the curve and the standard error
#   return(data.frame("n" = points_grid + 1, "pred" = pred, "se" = g_prime * se_fit))
# }




