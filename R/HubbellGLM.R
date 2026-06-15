#' Hubbell regression
#' This function mimics the GLM structure. We do so to avoid conflicts with the
#' current glm implementation (if it gets updated)
#'
#' @param formula Same as \code{stats::glm}
#' @param control Same as \code{stats::glm}
#' @param family a member of the hubbell family for \code{glm.fit.hubbell}.
#' @param data Same as \code{stats::glm}
#' @param weights Same as \code{stats::glm}
#' @param subset Same as \code{stats::glm}
#' @param na.action Same as \code{stats::glm}
#' @param start Same as \code{stats::glm}
#' @param etastart Same as \code{stats::glm}
#' @param mustart Same as \code{stats::glm}
#' @param offset Same as \code{stats::glm}
#' @param model Same as \code{stats::glm}
#' @param x Same as \code{stats::glm}
#' @param y Same as \code{stats::glm}
#' @param singular.ok Same as \code{stats::glm}
#' @param contrasts Same as \code{stats::glm}
#' @param ... Same as \code{stats::glm}
#'
#' @return An object of class \code{c('HubbellGLM', 'lm')}, which inherits from the class "lm"
#' @importFrom stats .getXlevels BIC add.scope coef deviance drop.scope extractAIC
#' @importFrom stats factor.scope fitted formula glm.fit is.empty.model
#' @importFrom stats model.extract model.frame model.matrix model.offset
#' @importFrom stats model.response model.weights na.pass napredict naresid
#' @importFrom stats nlminb nobs pnorm predict predict.lm pt residuals
#' @importFrom stats terms uniroot update update.formula vcov
#' @importFrom utils flush.console
#' @importFrom stats weights
#' @export
HubbellGLM <- function(formula,
                       family = hubbell(sigma = 0),
                       data,
                       weights,
                       subset,
                       na.action,
                       start = NULL,
                       etastart,
                       mustart,
                       offset,
                       control = list(...),
                       model = TRUE,
                       x = FALSE,
                       y = TRUE,
                       singular.ok = TRUE,
                       contrasts = NULL, ...) {
  # Force the method to be the glm for hubbell
  method <- "glm.fit.hubbell"

  # Here, the default glm code
  cal <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c(
    "formula", "data", "subset", "weights", "na.action",
    "etastart", "mustart", "offset"
  ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame")) {
    return(mf)
  }
  if (!is.character(method) && !is.function(method)) {
    stop("invalid 'method' argument")
  }
  if (identical(method, "glm.fit")) {
    control <- do.call("glm.control", control)
  }
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) {
      names(Y) <- nm
    }
  }
  X <- if (!is.empty.model(mt)) {
    model.matrix(mt, mf, contrasts)
  } else {
    matrix(, NROW(Y), 0L)
  }
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(weights) && any(weights < 0)) {
    stop("negative weights not allowed")
  }
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) {
      stop(gettextf(
        "number of offsets is %d should equal %d (number of observations)",
        length(offset), NROW(Y)
      ), domain = NA)
    }
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method,
    x = X, y = Y, weights = weights, start = start, etastart = etastart,
    mustart = mustart, offset = offset, family = family,
    control = control, intercept = attr(mt, "intercept") >
      0L, singular.ok = singular.ok
  ))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method,
      x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values,
      weights = weights, offset = offset, family = family,
      control = control, intercept = TRUE
    ))
    if (!fit2$converged) {
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    }
    fit$null.deviance <- fit2$deviance
  }
  if (model) {
    fit$model <- mf
  }
  fit$na.action <- attr(mf, "na.action")
  if (x) {
    fit$x <- X
  }
  if (!y) {
    fit$y <- NULL
  }

  # Return an object of class
  structure(c(fit, list(
    call = cal, formula = formula, terms = mt,
    data = data, offset = offset, control = control, method = method,
    contrasts = attr(X, "contrasts"),
    xlevels = .getXlevels(mt, mf)
  )), class = c(fit$class, c("HubbellGLM", "glm", "lm")))
}


#' @export
family.HubbellGLM <- function(object, ...) {
  object$family
}

#' @export
logLik.HubbellGLM <- function(object, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  fam <- family.HubbellGLM(object)$family
  dispersion <- family.HubbellGLM(object)$dispersion
  p <- object$rank
  val <- p - object$aic / 2
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}

#' @export
weights.HubbellGLM <- function(object, type = c("prior", "working"), ...)
{
  type <- match.arg(type)
  res <- if(type == "prior") object$prior.weights else object$weights
  if(is.null(object$na.action)) res
  else naresid(object$na.action, res)
}

#' @export
formula.HubbellGLM <- function(x, ...)
{
  form <- x$formula
  if( !is.null(form) ) {
    form <- formula(x$terms) # has . expanded
    environment(form) <- environment(x$formula)
    form
  } else formula(x$terms)
}


#' Estimate the sigma parameter for a HubbellGLM model via BIC minimisation
#'
#' @param formula A formula specifying the model (same syntax as
#'   \code{\link{HubbellGLM}}).
#' @param data A \code{data.frame} containing the variables in \code{formula}.
#' @param startpoint Numeric; starting value for the \code{nlminb} optimisation
#'   of \code{sigma}. Defaults to 0.
#' @param verbose Logical; if \code{TRUE} (default), prints the sigma value
#'   evaluated at each iteration.
#'
#' @return A scalar: the value of \code{sigma} that minimises the BIC of the
#'   fitted model.
#'
#' @export
estimate_sigma <- function(formula, data, startpoint = 0, verbose = TRUE){
  best_sigma <- nlminb(start = startpoint,
                       objective = function(x){
                         if (verbose) cat("Evaluating sigma = ", x, "\n")
                         fit <- HubbellGLM(formula,
                                           family = hubbell(link = "polyseries", sigma = x),
                                           data = data)
                         BIC(fit)
                       }, lower = -Inf, upper = 0.9, control = list(rel.tol = 1e-5))$par
  return(best_sigma)
}



# anova.HubbellGLM <- function (object, ..., dispersion = NULL, test = NULL)
# {
#
#   doscore <- !is.null(test) && test == "Rao"
#   varlist <- attr(object$terms, "variables")
#   x <- if (n <- match("x", names(object), 0L)){
#     object[[n]]
#   } else {
#     model.matrix(object)
#   }
#
#   varseq <- attr(x, "assign")
#   nvars <- max(0, varseq)
#   resdev <- resdf <- NULL
#   if (doscore) {
#     score <- numeric(nvars)
#     method <- object$method
#     y <- cbind(object$size, object$y)
#     fit <- eval(call(if (is.function(method)) "method" else method,
#                      x = x[, varseq == 0, drop = FALSE], y = y,
#                      weights = object$prior.weights,
#                      start = object$start, offset = object$offset, family = object$family,
#                      control = object$control))
#     r <- fit$residuals
#     w <- fit$weights
#     icpt <- attr(object$terms, "intercept")
#   }
#   if (nvars > 1 || doscore) {
#     method <- object$method
#     y <- cbind(object$size, object$y)
#     for (i in seq_len(max(nvars - 1L, 0))) {
#       fit <- eval(call(if (is.function(method)) "method" else method,
#                        x = x[, varseq <= i, drop = FALSE], y = y, weights = object$prior.weights,
#                        start = object$start, offset = object$offset,
#                        family = object$family, control = object$control))
#       if (doscore) {
#         zz <- eval(call(if (is.function(method)) "method" else method,
#                         x = x[, varseq <= i, drop = FALSE], y = r,
#                         weights = w, intercept = icpt))
#         score[i] <- zz$null.deviance - zz$deviance
#         r <- fit$residuals
#         w <- fit$weights
#       }
#       resdev <- c(resdev, fit$deviance)
#       resdf <- c(resdf, fit$df.residual)
#     }
#     if (doscore) {
#       zz <- eval(call(if (is.function(method)) "method" else method,
#                       x = x, y = r, weights = w, intercept = icpt))
#       score[nvars] <- zz$null.deviance - zz$deviance
#     }
#   }
#   resdf <- c(object$df.null, resdf, object$df.residual)
#   resdev <- c(object$null.deviance, resdev, object$deviance)
#   table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))),
#                       resdf, resdev)
#   tl <- attr(object$terms, "term.labels")
#   if (length(tl) == 0L)
#     table <- table[1, , drop = FALSE]
#   dimnames(table) <- list(c("NULL", tl), c("Df", "Deviance",
#                                            "Resid. Df", "Resid. Dev"))
#   if (doscore)
#     table <- cbind(table, Rao = c(NA, score))
#   title <- paste0("Analysis of Deviance Table", "\n\nModel: ",
#                   object$family$family, ", link: ", object$family$link,
#                   "\n\nResponse: ", as.character(varlist[-1L])[1L], "\n\nTerms added sequentially (first to last)\n\n")
#   df.dispersion <- Inf
#   if (is.null(dispersion)) {
#     dispersion <- summary(object)$dispersion
#     if (!is.null(object$family$dispersion) && is.na(object$family$dispersion)) {
#       df.dispersion <- object$df.residual
#     }
#   }
#   if (!is.null(test)) {
#     if (test == "F" && df.dispersion == Inf) {
#       fam <- object$family$family
#       if (fam == "binomial" || fam == "poisson")
#         warning(gettextf("using F test with a '%s' family is inappropriate",
#                          fam), domain = NA)
#       else warning("using F test with a fixed dispersion is inappropriate")
#     }
#     table <- stat.anova(table = table, test = test, scale = dispersion,
#                         df.scale = df.dispersion, n = NROW(x))
#   }
#   structure(table, heading = title, class = c("anova", "data.frame"))
# }




