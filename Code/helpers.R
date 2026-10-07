# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Shared helper functions for the modified variant surveillance pipeline
#
# This file holds the function definitions used by both
# variant_surveillance_modified.R and the unit tests, so that both point at
# the same code. It is derived from the original
# variant_surveillance_system.R (Paul, 2021-05-11 / 2022-01-18).
#
# MODIFICATION (relative to the original svymultinom):
#   The survey-design variance correction has been removed. svymultinom no
#   longer builds the score matrix or calls survey::svyrecvar, and it no
#   longer attaches `svyvcov` to the returned mlm object. As a result,
#   se.multinom automatically falls back to the model's own Hessian-based
#   covariance (solve(Hessian)) when computing standard errors. The point
#   estimates are unchanged (they come from the weighted nnet::multinom fit,
#   which never used the survey design). Only the prediction-interval
#   variance changes: it now reflects the weighted multinomial likelihood
#   treating observations as independent given the weights, rather than the
#   clustered/stratified survey design.
#
#   se.multinom and svyCI are unchanged from the original; they are included
#   here only so this file is self-contained.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

`%notin%` <- Negate(`%in%`)


# ---------------------------------------------------------------------------
# svymultinom (MODIFIED: Hessian-based variance, no survey-design correction)
# ---------------------------------------------------------------------------
svymultinom = function(mod.dat,
                       mysvy = NULL,
                       fmla = formula("as.numeric(as.factor(K_US)) ~ model_week"),
                       model_vars) {
  # Arguments:
  #  ~  mod.dat: source data (data.table)
  #  ~  mysvy:   (DEPRECATED) survey design object. No longer used for the
  #              variance correction; retained for signature compatibility.
  #              If supplied, its weights are used to aggregate the modeling
  #              data; if NULL, mod.dat$weights is used directly.
  #  ~  fmla:    multinomial model formula
  #  ~  model_vars: vector of variants included in the model (same order as K_US)
  #
  # Output: list with:
  #  mlm:       nnet multinomial model object (retains $Hessian; does NOT carry
  #             a $svyvcov, so se.multinom uses the Hessian-based covariance)
  #  estimates: coefficient estimates (named numeric vector)
  #  invinf:    variance-covariance matrix for coefficients = solve(-Hessian)
  #  SE:        SE of coefficients from the Hessian-based covariance
  #  variants:  modvars lookup table

  # get the variant names in the model
  modvars <- data.frame(
    'K_US' = 1:(length(model_vars)+1),
    'Variant' = c(names(model_vars), 'Other')
  )
  # get counts & weights by variant to help ID variants that might be causing issues
  moddatvars <- mod.dat[,
                        .(.N,
                          sum(weights)),
                        by = 'K_US']
  names(moddatvars) = c('K_US', 'N', 'Weight')
  modvars <- merge(x = modvars,
                   y = moddatvars,
                   all.x = TRUE,
                   by = 'K_US')

  # weights used to aggregate the modeling data. Prefer the survey design's
  # weights if a design is supplied (for backward compatibility); otherwise
  # use the weights column on mod.dat directly.
  agg_weights <- if (!is.null(mysvy)) weights(mysvy) else mod.dat$weights

  # aggregate data before fitting the multinomial model to improve run time
  fmla.vars = all.vars(fmla)
  mlm.dat = data.table::data.table(
    cbind(
      data.frame(mod.dat)[, fmla.vars],
      weight = agg_weights))[
        ,
        .(weight = sum(weight)), # aggregate "weight" column
        by = fmla.vars] # by the formula

  # Fit multinomial logistic regression
  # (without survey design, but with survey weights)
  multinom_geoid = nnet::multinom(formula = fmla,
                                  data    = mlm.dat,
                                  weights = weight,
                                  Hess    = TRUE,
                                  maxit   = 1000,
                                  trace   = FALSE)

  ## Format results
  # creates a list that contains the mlm object and the coefficient estimates
  rval = list(mlm = multinom_geoid,
              estimates = coefficients(multinom_geoid))

  # transforms the estimates object to be a list where each element is a vector
  #  of the coefficients for a given variant (us model: Intercept, model_week)
  rval$estimates = as.list(data.frame(t(rval$estimates)))

  # check if the Hessian is invertible; if not, return NA
  invinf <- tryCatch(
    {
      solve(-multinom_geoid$Hessian)
    },
    error = function(cond) {
      return(NA)
    }
  )

  # If the Hessian is NOT invertible, build a minimal output and warn.
  # se.multinom will still run but will assume SE is 0.
  if ( is.na(invinf[1]) ){

    # add empty items to match structure of a successful run
    rval <- append(rval,
                   list(
                     'invinf' = NULL,
                     'SE'     = NULL
                   ))

    # capture Hessian column names before removing the Hessian
    hess_colnames <- colnames(multinom_geoid$Hessian)
    hess_diag     <- diag(multinom_geoid$Hessian)

    # make the estimates a vector instead of a list
    rval$estimates = unlist(rval$estimates)

    # add prettier names to the estimates
    names(rval$estimates) <- paste0('b',
                                    sub(pattern = ':',
                                        replacement = '\\.',
                                        x = hess_colnames))

    # add a warning
    warning_message = 'Hessian is non-invertible. Nowcast estimates will not have prediction intervals. Check for variants with very few samples.'
    warning(warning_message)

    # print out counts by variant to help troubleshoot.
    # (USA-only: group by K_US only; the original grouped by HHS + K_US.)
    by_cols <- if ("HHS" %in% names(mod.dat)) c('HHS', 'K_US') else 'K_US'
    model_counts <- mod.dat[, sum(count), by = by_cols]
    model_counts$count = model_counts$V1
    model_counts$V1 = NULL
    model_counts = model_counts[order(model_counts$count, decreasing = FALSE), ]
    model_counts$Variant = modvars$Variant[model_counts$K_US]
    print('Here are counts by variant:')
    print(model_counts)

    # print out the highest and lowest values in the Hessian
    print('Also investigate very small or very large values in the Hessian')
    hess_headtail <- data.frame('element' = names(sort(hess_diag)),
                                'value' = sort(hess_diag))
    rownames(hess_headtail) <- 1:nrow(hess_headtail)
    print(hess_headtail[c(1:5, (nrow(hess_headtail)-4):nrow(hess_headtail)),])

    # NOTE: in the non-invertible case we leave the Hessian on the mlm object
    # set to NULL (as the original did) so se.multinom takes its zero-SE path.
    rval$mlm$Hessian = NULL

  } else {
    # MODIFIED: Hessian-based covariance, no survey-design correction.
    #
    # The original code built a score matrix from the multinomial gradient,
    # multiplied by invinf, and passed the result to survey::svyrecvar to get
    # a design-adjusted "sandwich" covariance, then attached it as
    # mlm$svyvcov. All of that is removed. We keep the inverse-information
    # (Hessian-based) covariance and report SEs from it. We deliberately do
    # NOT attach svyvcov, so se.multinom uses solve(mlm$Hessian) instead.

    # variance-covariance matrix for coefficients (not design-adjusted)
    rval$invinf = invinf

    # SE from the diagonal of the inverse-information matrix
    rval$SE = sqrt(diag(invinf))

    # convert the "estimates" object from a list to a vector
    rval$estimates = unlist(rval$estimates)

    # name the estimates to match the Hessian's coefficient names
    names(rval$estimates) <- paste0('b',
                                    sub(pattern = ':',
                                        replacement = '\\.',
                                        x = colnames(multinom_geoid$Hessian)))
    names(rval$SE) <- names(rval$estimates)

    # NOTE: mlm$Hessian is retained (not nulled), so se.multinom's Hessian
    # branch can compute solve(mlm$Hessian). No mlm$svyvcov is attached.
  }

  # return the rval object
  rval$variants = modvars
  return(rval)
}


# ---------------------------------------------------------------------------
# se.multinom (UNCHANGED from the original)
# ---------------------------------------------------------------------------
# With the modified svymultinom above, the returned mlm has a $Hessian but no
# $svyvcov, so the "Hessian" branch below is taken automatically. The sign
# convention (solve(Hessian), not solve(-Hessian)) is correct as documented.
se.multinom = function(mlm,
                       newdata_1row,
                       composite_variant = NULL) {
  # get the model coefficients
  cf = coefficients(mlm)

  # get the variance-covariance matrix
  if ("svyvcov" %in% names(mlm)) {
    vc = mlm$svyvcov
  } else {
    if ("Hessian" %in% names(mlm)) {
      vc = solve(mlm$Hessian)
      # solve(Hessian) (not -Hessian) because nnet/nlm minimizes the negative
      # log-likelihood, so the Hessian is already negated relative to the
      # log-likelihood Hessian.
    } else {
      # no covariance available -> zero variance
      vc = matrix(data = 0,
                  nrow = length(cf),
                  ncol = length(cf))
    }
  }

  # model matrix from the model fit & the new data
  mm = model.matrix(as.formula(paste("~", as.character(mlm$terms)[3])),
                    newdata_1row,
                    xlev = mlm$xlevels)

  # covariate names for each variant
  mnames = outer(X = mlm$lev[-1],
                 Y = colnames(mm),
                 FUN = paste, sep = ":")

  # matrix of covariate values
  cmat = matrix(data = 0,
                nrow = nrow(mnames),
                ncol = ncol(vc),
                dimnames = list(mlm$lev[-1],
                                as.vector(t(mnames))))
  for (rr in 1:nrow(cmat)) cmat[rr, mnames[rr,]] = c(mm)

  # linear predictor values
  y_i = c(0, coefficients(mlm) %*% c(mm))

  # variance-covariance matrix of the linear predictor
  vc.y_i = rbind(0, cbind(0, cmat %*% vc %*% t(cmat)))

  # predicted proportions
  p_i = exp(y_i) / sum(exp(y_i))

  # Taylor series based variance: dp_i/dy_j = delta_ij p_i - p_i * p_j
  dp_dy = diag(p_i) - outer(p_i, p_i, `*`)

  # variance/covariance of the predicted proportions
  p.vcov = dp_dy %*% vc.y_i %*% dp_dy
  se.p_i = as.vector(sqrt(diag(p.vcov)))

  # composite variants
  if (!is.null(composite_variant)) {
    composite_variant = list(
      matrix = composite_variant,
      p_i    = as.vector(composite_variant %*% p_i),
      se.p_i = as.vector(sqrt(diag(composite_variant %*% p.vcov %*% t(composite_variant))))
    )
  }

  # reshape vc to pull out the coefficient of time
  dim(vc) = rep(rev(dim(cf)), 2)

  return(
    list(p_i = p_i,
         se.p_i = se.p_i,
         b_i = c(0, cf[, 2]),
         se.b_i = c(0, sqrt(diag(vc[2,,2,]))),
         composite_variant = composite_variant)
  )
}


# ---------------------------------------------------------------------------
# svyCI (UNCHANGED from the original)
# ---------------------------------------------------------------------------
svyCI = function(p, s, ...) {
  if (s == 0) {
    return(c(0, 0))
  } else if (p == 0) {
    return(c(0, 0))
  } else if (p == 1) {
    return(c(1, 1))
  } else {
    n = p * (1 - p) / s^2
    out = prop.test(x = n * p,
                    n = n,
                    ...)$conf.int
    return(c(out[1], out[2]))
  }
}


# ---------------------------------------------------------------------------
# Lineage aggregation helpers: np / nearest_parent (UNCHANGED from original)
# ---------------------------------------------------------------------------
np = function(x, y, no_match = 'Other') {
  y <- y[ nchar(y) <= nchar(x) ]
  y. <- y
  y.[ nchar(y) < nchar(x) ] <- paste0(y[ nchar(y) < nchar(x) ], '.')
  sx = strsplit(x, "", fixed = TRUE)
  sy = strsplit(y., "", fixed = TRUE)
  match_array <- array(
    data = mapply(function(X, Y) {
      ylen <- seq_len(length(Y))
      wh <- (X[ylen] == Y[ylen])
      if (all(wh)) return(1) else (which.min(wh) - 1) / length(ylen)
    },
    rep(sx, each = length(sy)),
    sy),
    dim = c(length(x), length(y)),
    dimnames = list(x, y)
  )
  match100 <- colnames(match_array)[ match_array[1,] == 1 ]
  if (length(match100) > 0) {
    return(match100[ nchar(match100) == max(nchar(match100)) ])
  } else {
    return(no_match)
  }
}

nearest_parent <- function(x, y, no_match = 'Other') {
  sapply(x, function(x_i) np(x_i, y))
}
