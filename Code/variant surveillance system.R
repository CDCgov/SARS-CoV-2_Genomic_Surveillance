# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Code accompanying MMWR on Genomic Surveillance 
# Date: 2021-05-11; updated 2022-01-18
# Author: Prabasaj Paul
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Part 1: Define Functions ----------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Functions that do all the heavy lifting -------------------------------------

# Helper wrapper for survey::svyciprop
# this returns confidence intervals based on survey design only
# (no multinomial model for temporal smoothing)
myciprop = function(voc,
                    geoid,
                    svy,
                    str   = TRUE,
                    range = FALSE,
                    mut   = FALSE,
                    level = 0.95,
                    ci.type = 'KG',
                    ...){
  # Arguments:
  #  ~ voc:   character string of the variants of concern to analyze
  #           (if multiple listed, estimated proportion will be for the group of
  #           variants)
  #  ~ geoid: character string/numeric of the geographic resolution to analyze
  #           (options: 2 letter state code, 1:10 for HHS region, or "USA")
  #  ~ svy:   survey design object
  #  ~ str:   logical argument indicating if output should be formatted in single
  #           character object or as vector/dataframe
  #  ~ range: logical argument indicating whether to present CI range vs CI bounds
  #  ~ mut:   logical argument indicating whether to analyze mutation profiles or
  #           lineages
  #  ~ level: numeric specifying confidence level (1-alpha)
  #  ~ ...:   when method = "asin" or "mean", these arguments are passed on to
  #           svyciprop (which passes them to "confint.svystat")
  
  # Output: Confidence Intervals in 1 of several formats:
  #  1) if str == FALSE: named numeric vector with 7 values:
  #                      estimate: estimated proportion
  #                      lcl: lower confidence limit
  #                      ucl: upper confidence limit
  #                      DF: degrees of freedom
  #                      n.eff.capped: number of effective samples (capped to actual number sampled)
  #                      CVmean: coefficient of variation of the mean (estimated proportion)
  #                      deffmean: design effect of the mean (estimated proportion)
  #  2) if str == TRUE:
  #         if range == TRUE:
  #               string range: "ll.l-hh.h"
  #         if range == FALSE:
  #               string: "pp.p (ll.l-hh.h) DF=__, Eff.sampsize=__, CVprop=__, DEff=__"
  
  # Create a binary indicator variable to ID samples to analyze
  # (based on either S_MUT or VARIANT)
  if (mut) {
    VOC = grepl(
      pattern = paste0("^(?=.*\\",
                       paste(voc, collapse="\\b)(?=.*\\"),"\\b)"),
      x = svy$variables$S_MUT,
      perl = T
    )
  } else {
    VOC = (svy$variables$VARIANT %in% voc)
  }
  
  # extract the relevant survey design
  srv_design = subset(update(svy, VOC = VOC),
                      (STUSAB %in% geoid) |
                        (HHS %in% geoid) |
                        (geoid == "USA"))
  
  # calculate the confidence interval
  # (of all specified "voc" grouped together and all non-"voc" grouped together)
  if (ci.type == "xlogit"){
    # (note: this will throw a lot of warnings b/c of strata that have single
    # PSU's, so I'll suppress warnings for just this function)
    res = suppressWarnings(
      survey::svyciprop(~VOC,
                        design = srv_design,
                        method = "xlogit",
                        level = level,
                        ...)
    )
  } else {
    # use a modified version of svyciprop to limit effective sample size so that
    # it's never > observed sampled size.
    res = suppressWarnings(
      svycipropkg(~VOC,
                  design = srv_design,
                  level = level,
                  ...)
    )
  }
  
  # svyciprop & svycipropkg return a single number (the proportion estimate)
  # with the confidence interval in the attributes.
  # Add the confidence interval directly to the results
  res = c('estimate' = unname(res[1]),
          'lcl' = survey:::confint.svyciprop(res)[1],
          'ucl' = survey:::confint.svyciprop(res)[2])
  
  # add on the degrees of freedom
  res = c(res,
          'DF' = survey::degf(design = srv_design))
  
  
  #add effective sample size to res output
  m <- eval(bquote( # eval and bquote are superfluous b/c there's no variable substitution going on...
    # suppress warnings caused by single PSU in stratum
    suppressWarnings(
      #extracts the mean and SE for each estimate
      survey::svymean(x = ~as.numeric(VOC),
                      design = srv_design)
    )))
  
  # coefficient of variation
  CVmean <- suppressWarnings(
    # suppress warnings caused by single PSU in stratum
    # then calculate CV of the mean
    survey::cv(object = survey::svymean(x = ~as.numeric(VOC),
                                        design = srv_design))
  )
  
  # calculate the design effect
  # "The design effect compares the variance of a mean or total to the variance
  #  from a study of the same size using simple random sampling without replacement"
  deffmean <- suppressWarnings(
    # suppress warnings caused by single PSU in stratum
    survey::deff(object = survey::svymean(x = ~as.numeric(VOC),
                                          design = srv_design,
                                          deff = T))
  )
  
  # n_effective (this comes from the "survey::svyciprop" function when using the
  #              beta method for calculating the CI)
  #              coef(m) = mean     vcov(m) = variance
  n.eff <- coef(m) * (1 - coef(m))/vcov(m)
  
  # define alpha based on CI level
  alpha <- 1-level
  
  # this is also from the "survey::svyciprop" function using method == 'beta'
  n.eff.capped <- min(
    n.eff * (qt(p = alpha/2,
                df = nrow(srv_design) - 1) /
               qt(p = alpha/2,
                  df = survey::degf(design = srv_design)))^2,
    nrow(srv_design)
  )
  
  # add more info into the results
  res = c(res,
          'n.eff.capped' = n.eff.capped,
          'CVmean'= CVmean,
          'deffmean' = unname(deffmean))
  
  # optionally convert results to a string
  if (str) {
    res = c(round( 100 * res[1:3], 1),
            res[4],
            res[5],
            res[6],
            res[7])
    
    # optionally format results as a range
    if (range) {
      res = paste0(res[2], "-", res[3])
    } else {
      res = paste0(res[1],
                   " (",res[2],"-",res[3],
                   ") DF=",res[4],
                   ", Eff.sampsize=",res[5],
                   ", CVprop=",res[6],
                   ",DEff=",res[7])
    }
  }
  
  # return the results
  res
}

#Functions for multinomial nowcast model
# this function fits a multinomial model to the data to get temporally smoothed
# proportion estimates while also accounting survey design.
# Note: "svymultinom" is used in conjunction with "se.multinom" and "svyCI"
#       functions.
#       - "svymultinom" fits the model and adjusts the variance-covariance
#          matrix (of the model coefficients) for the sampling design.
#       - "se.multinom" uses the model coefficients from "svymultinom" to calculate
#          the linear predictor values (and covariances), and then converts the
#          linear predictor values into proportions and calculates the covariances
#          of those proportions.
#       - "svyCI" uses the proportions and SE from "se.multinom" to calculate
#          confidence intervals on estimated proportions (using score intervals
#          for a Binomial distribution rather than the (less reliable) normal
#          approximation)
svymultinom = function(mod.dat,
                       mysvy,
                       fmla = formula("as.numeric(as.factor(K_US)) ~ model_week + as.factor(HHS)"),
                       model_vars) {
  # Arguments:
  #  ~  mod.dat: source data (data.table)
  #  ~  mysvy:   survey design object
  #  ~  fmla:    multinomial model formula
  #  ~  model_vars: vector of variants that are included in the model (in same order as K_US). 
  #                 This is only used for helping to troubleshoot model-fit issues. 
  
  # Output: list of 6 objects:
  #  mlm:       nnet multinomial model object (edited to include "svyvcov", which
  #             is the variance-covariance matrix after accounting for survey
  #             design (i.e. "sandwich" below))
  #  estimates: coefficient estimates from the mlm model (named numeric vector)
  #  scores:    gradient of the loglikelihood * survey weights * model matrix
  #             (numeric matrix with row for each observation and column for each coefficient)
  #  invinf:    variance-covariance matrix for coefficients of the mlm model (i.e.
  #             not accounting for survey design)
  #             (numeric matrix with row and column for each coefficient)
  #  sandwich:  variance-covariance matrix for coefficients after accounting for
  #             survey design
  #             (numeric matrix with row and column for each coefficient)
  #  SE:        SE of coefficients after accounting for survey design
  #             (named numeric vector)
  
  
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
  
  # aggregate data before fitting the multinomial model to improve run time
  fmla.vars = all.vars(fmla)
  mlm.dat = data.table::data.table(
    cbind(
      data.frame(mod.dat)[, fmla.vars],
      weight = weights(mysvy)))[
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
  
  ## Format results to fit into svymle-like workflow
  # get the number of variants being modeled
  # (should be length of model_vars plus 1 for others)
  num_var = length(unique(with(mod.dat,
                               eval(terms(fmla)[[2]]))))
  
  # generate the model formula without the response term
  fmla.no.response = formula(delete.response(terms(fmla)))
  
  # repeat the model formula w/o response term for each variant listed in model_vars
  formulas = rep(list(fmla.no.response),
                 num_var - 1)
  
  # Add response back to first formula to ensure inclusion as response later on
  formulas[[1]] = fmla
  
  # sets the names of the formulas to correspond to beta coefficients for each variant
  names(formulas) = paste0("b", 1:length(formulas) + 1)
  
  # creates a list that contains the mlm object and the coefficient estimates
  # (this will form the output of this function after other things are added on)
  rval = list(mlm = multinom_geoid,
              estimates = coefficients(multinom_geoid))
  
  # transforms the estimates object to be a list where each element is a vector
  #  of the coefficients for a given variant
  # (hhs model: Intercept, week, HHS 2:10;
  #   us model: Intercept, week)
  rval$estimates = as.list(data.frame(t(rval$estimates)))
  ## End multinomial regression
  
  # check if the Hessian is invertible
  #   if not, return NA
  invinf <- tryCatch( 
    { 
      solve(-multinom_geoid$Hessian) 
    }, 
    error = function(cond) { 
      return(NA) 
    } 
  ) 

  # If the Hessian is NOT invertible, just create the output.
  # (note: the se.multinom function will still run on the output of this function,
  #  but will assume the SE is 0)
  if ( is.na(invinf[1]) ){ 
    
    # add empty items to the results to match structure of a successful run
    # (setting to NULL inside of list() *does* create the named item)
    rval <- append(rval,
                   list(
                     'scores'   = NULL,
                     'invinf'   = NULL,
                     'sandwich' = NULL,
                     'SE'       = NULL
                   ))
    
    # remove the Hessian
    # (setting to NULL outside of list() *removes* item)
    rval$mlm$Hessian = NULL
    
    # make the estimates a vector instead of a list
    rval$estimates = unlist(rval$estimates)
    
    # add prettier names to the estimates (to match names that are produced below)
    names(rval$estimates) <- paste0('b',
                                    sub(pattern = ':',
                                        replacement = '\\.',
                                        x = colnames(multinom_geoid$Hessian)))
    
    # add a warning
    warning_message = 'Hessian is non-invertible. Nowcast estimates will not have prediction intervals. Check for geographic regions with very few samples of a variant.'
    warning(warning_message)
    
    # print out regions with the fewest observations
    # if there are problems with the Hessian, investigate counts by region
    model_counts <- mod.dat[,sum(count),by = c('HHS', 'K_US')]
    model_counts$count = model_counts$V1
    model_counts$V1 = NULL
    model_counts = model_counts[order(model_counts$count, decreasing = FALSE), ]
    model_counts$Variant = modvars$Variant[model_counts$K_US]
    print('Here are counts by region:')
    print(model_counts[,c('HHS', 'K_US', 'Variant', 'count')])
    
    # print out the highest and lowest values in the Hessian
    print('Also investigate very small or very large values in the Hessian')
    hess_headtail <- data.frame('element' = names(sort(diag(multinom_geoid$Hessian))),
                                'value' = sort(diag(multinom_geoid$Hessian)))
    rownames(hess_headtail) <- 1:nrow(hess_headtail)
    print(hess_headtail[c(1:5, (nrow(hess_headtail)-4):nrow(hess_headtail)),])
    
  } else {
    # if the Hessian is solvable, then adjust variance-covariance matrix using survey design
    
    ## Define likelihood
    
    # Loglikelihood:
    lmultinom = function(v, ...) {
      # vectorized loglikelihood function
      # Arguments: v  (positive integer vector) is position of response variable
      #                in ordered list of possible values (1=reference)
      #            ...  vectors of  linear predictors
      
      # create a matrix of linear predictors (with 0's for reference class)
      b = cbind(0, ...)
      b = matrix(b,
                 nrow = length(v))
      
      # get the predicted probability of the realized (actual) outcomes
      sapply(X = seq_along(v),
             FUN = function(rr) b[rr, v[rr]])  - log(rowSums(exp(b)))
    }
    
    # Gradient of loglikelihood:
    gmultinom = function(v, ...) {
      # vectorized partial derivatives of lmultinom(v, ...) with respect to
      #  linear predictors at ...
      # vectorized loglikelihood function
      # Arguments: v  (positive integer vector) is position of response variable
      #                in ordered list of possible values (1=reference)
      #            ...  vectors of linear predictors
      
      # create a matrix of linear predictors (with 0's for reference class)
      #  column for each observation
      #  row for each possible outcome class
      b = cbind(0, ...)
      b = matrix(b,
                 nrow=length(v))
      
      # convert to probabilities by normalizing linear predicted values by the
      # total value across all possible outcome classes
      p = exp(b) / rowSums(exp(b))
      
      # create a matrix of the same dimensions
      # to hold the observed class of each observation
      delta_ij = p * 0
      
      # set some values of the matrix to 1
      # (for each observation (i.e. row) assign a value of 1 to the column
      #  representing its clade/rank order)
      delta_ij[1:length(v) + (v-1) * length(v)] = 1
      
      # Column n is partial derivative with respect to ...[n]:
      # (realized outcome minus the predicted probability of belonging to that class)
      (delta_ij - p)[, -1]
    }
    ## End likelihood definitions
    
    ## Build dataframe to pass to svyrecvar
    # Adapted from code in svymle
    # modify formula names
    nms = c("", names(formulas))
    
    # logical that identifies which formula has a response variable
    has.response = sapply(formulas, length) == 3
    
    # creates variable equal to regression formula that has response variable
    ff = formulas[[which(has.response)]]
    
    #sets predictor terms to 1 so: as.numeric(as.factor(K_US)) ~ 1
    ff[[3]] = 1
    
    #I think this just gets you what the rankings are for the variants from mod.dat
    # (i.e. the response data for the regression model)
    y = eval.parent(model.frame(formula   = ff,
                                data      = mod.dat,
                                na.action = na.pass))
    
    # sets the first formula to a formula without the response term
    formulas[[which(has.response)]] = formula(delete.response(terms(formulas[[which(has.response)]])))
    
    # vector with the names of the regression predictors
    vnms = unique(do.call(c, lapply(formulas, all.vars)))
    
    # converts the vnms to a formula object w/o the response term
    uformula = make.formula(vnms)
    
    # gets the values for the model predictors from mod.dat
    mf = model.frame(formula   = uformula,
                     data      = mod.dat,
                     na.action = na.pass)
    
    # adds a column for the variant rankings to mf with the column header as
    #  as.numeric(as.factor(K_US))
    mf = cbind(`(Response)` = y,
               mf)
    
    # ensure that there aren't any duplicated columns
    #  (the drop = FALSE ensures that mf remains a dataframe even if only 1 column
    #   is returned)
    mf = mf[, !duplicated(colnames(mf)),
            drop = FALSE]
    
    # get the svy weights from the survey design object
    weights = weights(mysvy)
    
    # defines Y (response variable) as variant group (rank abundance)
    Y = mf[, 1]
    
    # for each formula listed in "formulas", generate the underlying data for
    # those parameters
    # (model.frame object is required for the model.matrix function)
    mmFrame = lapply(X    = formulas,
                     FUN  = model.frame,
                     data = mf) # argument passed to "model.frame"
    
    # creates model design objects for each rank class
    mm = mapply(FUN      = model.matrix,
                object   = formulas,
                data     = mmFrame,
                SIMPLIFY = FALSE)
    
    # get the cumulative number of columns across all dataframes listed in mm
    np = c(0,
           cumsum(sapply(X   = mm,
                         FUN = NCOL)))
    
    # get the column names from all the dataframes in mm
    parnms = lapply(mm, colnames)
    
    # format column names in each list element to correspond to the coefficient
    for (i in 1:length(parnms)){
      parnms[[i]] = paste(nms[i + 1],
                          parnms[[i]],
                          sep = ".")
    }
    
    # convert the parameter names from a list to a vector
    parnms = unlist(parnms)
    
    # get the estimated coefficients from multinom as a vector
    theta = unlist(rval$estimates)
    
    # creates empty list the same length as the number of variants plus "other"
    args = vector("list", length(nms))
    
    # sets first list element to the variant rank data (i.e. response variable)
    args[[1]] = Y
    
    # assigns the names of the args list to nms
    names(args) = nms
    
    # get the linear predictor for each response variable class (rank abundance)
    # by multiplying the covariate values in "mm" by the coefficient values, "theta"
    for (i in 2:length(nms)){
      args[[i]] = drop(mm[[i - 1]] %*% theta[(np[i - 1] + 1):np[i]])
    }
    
    # this takes the args list and uses the gmultinom function above to estimate
    # the gradient of the log likelihood
    deta = matrix(data = do.call(what = "gmultinom",
                                 args = args),
                  ncol = length(args) - 1)
    
    # creates an empty list element called scores in the rval list
    rval <- append(rval,
                   list('scores' = NULL))
    
    # creates numerical vector the same length as the names of formulas
    reorder = na.omit(match(x = names(formulas),
                            table = nms[-1]))
    
    #for each variant this multiplies the gradient of the loglikelihood by the survey weights
    for (i in reorder){
      rval$scores = cbind(rval$scores,
                          deta[, i] * weights * mm[[i]])
    }
    
    # invert negative Hessian to get an estimate of variance-covariance matrix
    rval$invinf = invinf 
    
    # assign names to the rval$invinf
    dimnames(rval$invinf) = list(parnms, parnms)
    
    # Use matrix multiplication to multiply the "scores" (i.e., the gradient
    # loglikelihood * survey weights * model matrix) by the variance-covariance
    # matrix of parameter estimates
    db = rval$scores %*% rval$invinf
    
    # Everything up until now has been formatting steps to get the estimates/data
    # formatted in way that's compatible with the svyrecvar function
    
    # Computes the variance of a total under multistage sampling, using a recursive
    # descent algorithm. The object that's returned ("sandwich") is the covariance matrix
    rval$sandwich = survey::svyrecvar(x          = db,
                                      clusters   = mysvy$cluster,
                                      stratas    = mysvy$strata,
                                      fpcs       = mysvy$fpc,
                                      postStrata = mysvy$postStrata)
    
    # estimate the SE as the square root of the diagonal elements of the
    # variance-covariance matrix (diagonal = variance)
    rval$SE = sqrt(diag(rval$sandwich))
    
    # convert the "estimates" object from a list to a vector
    rval$estimates = unlist(rval$estimates)
    
    # assign the names of "estimates" in rval to be the same as the names of "SE"
    names(rval$estimates) = names(rval$SE)
    
    # adds the variance-covariance matrix to the mlm object (mlm being the results
    # object you get from solving the multinomial model)
    rval$mlm$svyvcov = rval$sandwich
  }
  
  # return the rval object
  rval$variants = modvars
  return(rval)
}


# function to extract the predictions and SE of predictions from the multinomial
# nowcast model for a particular week and location
# (and optionally grouped clades/variants)
# Note: This function estimates proportions (predicted values) and their standard
#       errors using the results of "svymultinom"; the "svymultinom" function
#       fits the model and estimates the SE of model coefficients.
se.multinom = function(mlm,
                       newdata_1row,
                       composite_variant = NULL) {
  # Arguments:
  #  ~  mlm:   model output, with Hessian;
  #  ~  newdata_1row:  1-row dataframe consistent with predictors in formula
  #                    typically includes focal week (model_week) & geographic region
  #  ~  composite_variant: (if not NULL) is a matrix with one column per model
  #                        lineage and one row per composite variant. Matrix
  #                        element of 1 marks each model lineage (column)
  #                        for aggregation into a composite variant (row).
  #                        Example: matrix(c(1, 0, 0, 1, 0, 0), nrow = 1)
  #                        will combine the first and fourth lineages in the 
  #                        model into a single estimate, which will be appended
  #                        at end of p_i and se.p_i
  
  # If mlm is without Hessian, all se's are set to zero
  # mlm not geographically stratified for geoid="USA": ~ (week - current_week)
  # mlm, for HHS Regions: ~ week + HHS (1 as reference level)
  
  # Output: list of 5 items:
  #  ~  p_i:     predicted values for each clade/variant (numeric vector)
  #  ~  se.p_i:  SE of predicted values (numeric vector)
  #  ~  b_i:     coefficient (beta) values (numeric vector)
  #              NOTE! b_i is only interpretable as growth rates if first
  #                    predictor term is time!
  #  ~  se.b_i:  SE of coefficients (numeric vector)
  #  ~  composite_variant: list with 3 element:
  #                        matrix = composite_variant (input matrix)
  #                        p_i    = predicted proportions for composite variants
  #                        se.p_i = SE of the predicted proportions
  
  # get the model coefficients
  cf = coefficients(mlm)
  
  # get the variance-covariance matrix
  if ("svyvcov" %in% names(mlm)) {
    vc = mlm$svyvcov
  } else {
    if ("Hessian" %in% names(mlm)) {
      vc = solve(mlm$Hessian)
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      # solve(hessian) instead of solve(-hessian) here because
      #   svymle(), through a call to nlm(), minimizes a function. The functions
      #   lmultinom() and gmultinom() are, therefore, negative log-likelihood
      #   and it's gradient, respectively. The Hessian is, therefore, negative
      #   of what you'd expect.
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    } else {
      # If there's no variance-covariance matrix from the svymultinom function
      # and there's no Hessian, just set the variance-covariance to 0
      vc = matrix(data = 0,
                  nrow = length(cf),
                  ncol = length(cf))
    }
  }
  
  # get the model matrix from the model fit & the new data
  mm = model.matrix(as.formula(paste("~",as.character(mlm$terms)[3])),
                    newdata_1row,
                    xlev = mlm$xlevels)
  
  # get the covariate names for each variant
  mnames = outer(X = mlm$lev[-1],
                 Y = colnames(mm),
                 FUN = paste, sep=":")
  
  # matrix of covariate values (first create empty matrix)
  cmat = matrix(data = 0,
                nrow = nrow(mnames),
                ncol = ncol(vc),
                dimnames = list(mlm$lev[-1],
                                as.vector(t(mnames))))
  
  # fill in covariate values into the cmat
  for (rr in 1:nrow(cmat)) cmat[rr, mnames[rr,]] = c(mm)
  
  # linear predictor values
  y_i = c(0, coefficients(mlm) %*% c(mm))
  
  # get the variance-covariance matrix of the linear predictor
  # (after adjusting for survey design)
  vc.y_i = rbind(0, cbind(0, cmat %*% vc %*% t(cmat)))
  
  # calculate the predicted proportions
  p_i = exp(y_i)/sum(exp(y_i))
  
  # Taylor series based variance:
  # dp_i/dy_j = delta_ij p_i - p_i * p_j
  dp_dy = diag(p_i) - outer(p_i, p_i, `*`)
  
  # calculate variance/covariance of the predicted proportions
  p.vcov = dp_dy %*% vc.y_i %*% dp_dy
  # SE = sqrt of variance
  se.p_i = as.vector(sqrt(diag(p.vcov)))
  
  # if there are composite variants, calculate calculate their proportions & SE
  if (!is.null(composite_variant)) {
    composite_variant = list(
      matrix = composite_variant,
      p_i    = as.vector(composite_variant %*% p_i),
      se.p_i = as.vector(sqrt(diag(composite_variant %*% p.vcov %*% t(composite_variant))))
    )
  }
  
  # change the dimensions of the variance-covariance matrix to aid in
  # pulling out the coefficient of time.
  dim(vc) = rep(rev(dim(cf)), 2)

  # output
  return(
    list(p_i = p_i, # predicted values
         se.p_i = se.p_i, # SE of predicted values
         b_i = c(0, cf[, 2]), # coefficient value for week
         se.b_i = c(0, sqrt(diag(vc[2,,2,]))), # SE of the coefficient
         composite_variant = composite_variant) # predicted probabilities (and SE) for the composite variants
  )
}

# Function to get binomial confidence interval based on the point estimated
# proportion and associated SE
# (based on output from svymultinom & se.multinom functions)
svyCI = function(p, s, ...) {
  # p = point estimate
  # s = estimated standard error
  # ... optional arguments passed on to prop.test (e.g. conf.level = 0.95)
  
  # if se is 0, n will be Inf (possibly because of non-invertible Hessian); return CI of 0,0
  if (s == 0) {
    return(c(0, 0)) # return confidence interval of [0,0]
  } else if (p == 0) {
    # if p == 0 or p == 1, then n will be 0 & prop.test will throw an error
    return(c(0, 0)) # return confidence interval of [0,0]
  } else if (p == 1) {
    return(c(1, 1)) # return confidence interval of [1,1]
  } else {
    # calculate the sample size
    n = p*(1-p)/s^2
    
    # calculate the confidence interval using a proportion test
    out = prop.test(x = n * p,   # a vector of counts of successes
                    n = n,       # a vector of counts of trials
                    ...          # optional arguments (e.g. conf.level = 0.99)
    )$conf.int

    # return the lower and upper confidence interval limits
    return(c(out[1],out[2]))
  }
}


# svycipropkg function – Korn and Graubard confidence limits for proportion from complex survey data
# based on svyciprop function with "beta" method from survey package
# 'design' can be an age standardized object created using 'svystandardize'
# this implementation: Crescent Martin, NCHS
# last updated 6/3/2019
#NOTE: the svyciprop function from the survey package uses the calculated
#     effective sample size for the KG CI (which can be larger than the actual sample size)
#     The svycipropkg function caps the effective sample size at the actual sample size
#CAVEATE: the svycipropkg code hasn't gone through comprehensive testing.
svycipropkg <- function(formula, 
                        design, 
                        level = 0.95, 
                        df = survey::degf(design), 
                        ...) {
  
  # based on svyciprop function in the survey package with method = "beta"
  m <- eval(bquote(svymean(~as.numeric(.(formula[[2]])), design, ...)))
  rval <- coef(m)[1]
  attr(rval, "var") <- vcov(m)
  alpha <- 1 - level
  
  # if design is an age-standardized survey design object - custom code not in svyciprop function
  if (!is.null(design[['postStrata']])){
    # error checking
    if (!as.character(design[['call']][[1]]) == "svystandardize"){
      stop("svycipropkg design cannot be a subset of an age-standardized survey design object")
    }
    
    # create temporary data frame  
    design_tmp <- design[['variables']][which(design[['prob']] != Inf), ]
    # get the age-adjustment variable
    ageadjvar <- as.character(design[['call']]$by[[2]])
    # get the age-adjustment population weights
    population <- eval(design[['call']]$population)
    # calculate the weighted proportion and get the sample size for each age adjustment group
    p <- lapply(split(design_tmp, design_tmp[[ageadjvar]]),
                function(x) weighted.mean(x[[as.character(formula[[2]])]],
                                          x[[names(design[['allprob']])]] ) )
    n <- lapply(split(design_tmp, design_tmp[[ageadjvar]]), nrow )
    p <- unlist(p)
    n <- unlist(n)
    # calculate the variance of a simple random sample of the same size, for each age group
    age_var <- p*(1-p)/n
    # normalize the age-adjustment weights to total 1
    pop <- population/sum(population)
    # accumulative the SRS variance over age groups
    varsrs_adj <- sum(pop^2 * age_var)
    # design effect
    deff_adj <- ifelse(sum(p) == 0, 1, vcov(m)/varsrs_adj)
    # adjusted effective sample size
    n.eff <- ifelse(rval == 0, nrow(design_tmp),
                    min(nrow(design_tmp),
                        nrow(design_tmp)/deff_adj * (qt(alpha/2, nrow(design_tmp) - 1)/qt(alpha/2, degf(design)))^2))
    
  }
  else { # crude estimates (not age-adjusted)
    # effective sample size
    n.eff <- coef(m) * (1 - coef(m))/vcov(m)
    # modification to svyciprop: cap adjusted effective sample size at the actual sample size
    # Under the assumption that the true design effect is >1, though it may be estimated as <1 due to instability of the variance estimator 
    # This modification produces different estimated CIs than svyciprop IF the estimated design effect <1
    n.eff <- min( n.eff * (qt(alpha/2, nrow(design) - 1)/qt(alpha/2, degf(design)))^2, nrow(design))
  }
  
  ci <- c(qbeta(alpha/2, n.eff * rval, n.eff * (1 - rval) +  1), qbeta(1 - alpha/2, n.eff * rval + 1, n.eff * (1 - rval)))
  halfalpha <- (1 - level)/2
  names(ci) <- paste(round(c(halfalpha, (1 - halfalpha)) * 
                             100, 1), "%", sep = "")
  names(rval) <- paste(deparse(formula[[2]]), collapse = "")
  attr(rval, "ci") <- ci
  class(rval) <- "svyciprop"
  rval
}

# define "%notin%" function
`%notin%` <- Negate(`%in%`)


# define some functions for lineage aggregation
{
  # define a function to find a variant's nearest parent from a list of
  # options. It defines "parent" as the longest string that is a subset
  # of the child variant name. (So it only works on extended names, not on Pango abbreviations.)
  # If no match is found, it returns "Other"
  # takes 3 arguments:
  # arg x  = (string) variant name to look up
  # arg y  = (string) vector of voc (extended)
  # arg no_match = (string) what to return if there is no match
  # returns string with the parent's name
  np = function(x, y, no_match = 'Other') {

    # only look for parents in y, so each y MUST be shorter than x
    y <- y[ nchar(y) <= nchar(x) ]

    # add a "." to all y that are shorter than x (to make sure that the parent variant IS a parent)
    # e.g. without the extra ".", the parent of BA.5.22 might be identified as BA.5.2 (shorter & a perfect subset)
    y. <- y
    y.[ nchar(y) < nchar(x) ] <- paste0(y[ nchar(y) < nchar(x) ], '.')

    # split out each character in x & y
    sx = strsplit(x, "", fixed = TRUE)
    sy = strsplit(y., "", fixed = TRUE)

    # calculate an array of match proportions
    match_array <- array(
      # cycle over each string in X and Y
      data = mapply(function(X, Y) {
        # get the length of the shorter (parent) string
        #slen = seq_len(min(length(X), length(Y)))
        ylen <- seq_len(length(Y))

        # test if the first slen characters are the same in both strings
        # if they're not the same, calculate the proportion of characters (starting from the beginning) that are the same (i.e. how far you get through the string before the first mismatch)
        # wh <- (X[slen] == Y[slen])
        # if(all(wh)) return(1) else (which.min(wh) - 1) / length(slen)
        wh <- (X[ylen] == Y[ylen])
        if(all(wh)) return(1) else (which.min(wh) - 1) / length(ylen)
      },
      # values of X to pass to mapply function (just repeat x for each item in y)
      rep(sx, each = length(sy)),
      # values of Y to pass to mapply function
      sy),
      # dimensions for the array
      dim = c(length(x), length(y)),
      # names of the dimensions for the array
      dimnames = list(x, y)
    )
    # match_array has x strings on the x-axis and y strings on the
    # y-axis. The numbers in the array are the proportion of the 2 strings
    # that match (starting from beginning; so the proportion is how far through
    # the string you get before the first mismatch). A value of 1
    # means that the shorter string is a subset of the longer string.

    # id which matches are complete matches
    match100 <- colnames(match_array)[ match_array[1,] == 1 ]

    # if there are complete matches, get the parent
    if( length(match100) > 0 ) {

      # the longest complete match is the closest parent
      return(match100[ nchar(match100) == max(nchar(match100)) ])

      # # return the longest complete match (that is still shorter than x) as the nearest parent
      # match100_2 <- match100[ (nchar(match100) < nchar(x)) ]
      #
      # # if there are complete matches where
      # if( length(match100_2) > 0 ) return(match100_2[ nchar(match100_2) == max(nchar(match100_2)) ])
      # else return(no_match)
    } else {
      # if there are no complete matches, return "Other"
      return(no_match)
    }
  }
  # sapply(sort(unname(unlist(extra_voc_to_consider))), function(x) nearest_parent(x, vocxl))

  # convert function "np" to work on a vector of x (reuses the same y vector for each element in x)
  nearest_parent <- function(x, y, no_match = 'Other') {
    # this just uses "sapply" to look for the nearest parent of each x in y
    sapply( x, function(x_i) np(x_i,  y ) )
  }

  # # This is a shorter & easier-to-understand version of the "np" function
  # # (but it runs about 8 times slower than "np")
  # np2 <- function(x, y, no_match = 'Other'){
  #
  #    # y is a parent of x if x & y start with identical characters, then the next character in x is either "\\." or end of string "$"
  #    # all the parent variants in y
  #    all_pv <- y[unlist(lapply(y, function(p){
  #       grepl(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^', p)), '((\\.)|($))'),
  #            x = x)
  #    }))]
  #
  #    # if there are complete matches, get the parent
  #    if( length(all_pv) > 0 ) {
  #
  #       # only the max-length parent variant
  #       all_pv[ nchar(all_pv) == max(nchar(all_pv))]
  #
  #    } else {
  #       # if there are no complete matches, return "Other"
  #       return(no_match)
  #    }
  # }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Part 2: Data Prep & weight calculation --------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library(survey)     # package with survey design functions
library(nnet)       # package with multinomial regression for "Nowcast" model
library(data.table) # package for speeding up calculations

# set global options:
#  - tell the survey package how to handle surveys where only a single primary
#    sampling units has observations from a particular domain or subpopulation.
#    See more here:
#    https://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
#  - do not convert strings to factors (as was default in R < 4.0)
options(survey.adjust.domain.lonely = T,
        survey.lonely.psu = "average",
        stringsAsFactors = FALSE)

## Load sequence data with test tallies ----------------------------------------
dat <- readRDS(file = 'example_data.RDS')
# these data have already been filtered to ensure: 
# 1) all sequences are from human hosts
# 2) all sequences originated from the USA
# 3) all sequences have valid state abbreviations
# 4) there are no duplicates
# 5) all dates are valid (i.e. >= 2019-10-01 and <= (day of the analysis))
# 6) all sequences come from labs that provided at least 100 sequences total
# 7) all sequences are attributable to a lab (i.e. drop NA's)
# 8) all sequences have a valid variant name (ie. drop NA's and "None")

# create a variable for aggregating counts of different groupings of samples
dat$count = 1

# convert the data.frame to a data.table for faster calculation
dat = data.table::as.data.table(dat)



## Some general parameters -----------------------------------------------------

# starting date for defining weeks
week0day1 = as.Date("2020-01-05")

# Variants to be included in the analysis (not limited to VOCs; can be any Pango lineage)
voc = c("BA.1.1",
        "BA.2",
        "BA.2.12.1",
        "BA.2.75",
        "BA.2.75.2",
        "CH.1.1",
        "BF.7",
        "BF.11",
        'BA.4',
        'BA.4.6',
        'BA.5',
        'BA.5.2.6',
        "BQ.1",
        "BQ.1.1",
        "BN.1",
        "XBB",
        'XBB.1.5',
        'XBB.1.5.1',
        'XBB.1.9.1',
        'FD.2',
        'XBB.1.9.2',
        'XBB.1.16',
        'XBB.2.3',
        "B.1.617.2", # Delta
        "B.1.1.529") # Omicron

# fewer variants runs faster
# voc = voc[16:18]

## Aggregate sublineages -------------------------------------------------------
# This is an inefficient system based on abbreviated Pango lineage names, which necessitates a lot of code.
# All that it does is assign each Pango lineage to its nearest parent lineage that is included in "voc".
# make sure "VARIANT" is a character (rather than factor)
# (redundant with "stringsAsFactors = FALSE")
dat$VARIANT <- dat$lineage <- as.character(dat$VARIANT)

# create a table of abbreviated & expanded lineages
lut <- dat[, .(expanded_lineage = unique(expanded_lineage)), by = 'VARIANT']

# lineage_expanded for each voc
voc_expanded <- lut[ VARIANT %in% voc ][['expanded_lineage']]

# all variants in the data
unique_vars <- na.omit(unique(dat$VARIANT))
# lineage_expanded for each of unique_vars
unique_lineage_expanded <- unname(setNames(lut$expanded_lineage, lut$VARIANT)[ unique_vars ])

# create a lookup table for all the variants in the data
# (originally "lut" included all Pango lineages in the database)
# variant = abreviated variant name
# lineage_expanded
# parent_lineage_expanded = lineage_expanded that this variant will be aggregated into
# parent_variant = variant that this variant will be aggregated into
voc_lut <- data.frame(
'variant'                  = unique_vars,
'lineage_expanded'         = unique_lineage_expanded,
'parent_lineage_expanded'  = nearest_parent(unique_lineage_expanded, voc_expanded)
)

# add in the parent variant
voc_lut$parent_variant <- setNames(voc_lut$variant, voc_lut$lineage_expanded)[voc_lut$parent_lineage_expanded]
# update the row names
row.names(voc_lut) <- 1:nrow(voc_lut)


# use the voc_lut to get the aggregation variant ("parent_variant") for each variant
dat$VARIANT <- setNames(voc_lut$parent_variant, voc_lut$variant)[ dat$VARIANT ]
# see which variants were changed: data.table:::unique.data.table(dat[VARIANT != lineage, .(VARIANT, lineage)])

# create another column for the varients of interest
# this is only used to get (unweighted) counts of the sequences by lineage (used in all runs)
dat$VARIANT2 = as.character(dat$VARIANT)
# group all non-"variants of interest" together
dat[dat$VARIANT %notin% voc, "VARIANT2"] <- "Other"


# SGTF over-sampling was a temporary issue and has been removed



## Survey Weights --------------------------------------------------------------
# Estimation of the infection weights involves estimating the (unobserved) number
# of infections from test results. There is no reliable and precise method for this
# yet, so these weights are subject to considerable uncertainties.
# We use methods in:
#     https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm
#     https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009374
# Other strategies are also available; e.g.
#     https://covid19-projections.com/estimating-true-infections-revisited/


# UPDATED weights using NREVSS testing data
# proxy_infections    = estimated number of infections in a given state-week
# state_population    = state population
# POSITIVE.HHS.nrevss = number of positive tests reported to NREVSS in a given HHS region in a given week
# TOTAL.HHS.nrevss    = number of total    tests reported to NREVSS in a given HHS region in a given week
# POSITIVE.HHS        = number of positive tests reported to CLERS in a given HHS region in a given week
# population_reporting.HHS = total population of the states in a given HHS region that reported tests to CLERS in a given week
dat[
   ,
   "proxy_infections" := state_population * sqrt(POSITIVE.HHS.nrevss / TOTAL.HHS.nrevss) *
      sqrt( POSITIVE.HHS / population_reporting.HHS )]
# this formula is equivalent to calculating the regional infections as: sqrt(POSITIVE.HHS.nrevss / TOTAL.HHS.nrevss * population_reporting.HHS * POSITIVE.HHS) and then splitting it up among states in the region based on population (i.e. multiplying by): (state_population / population_reporting.HHS)

# calculate the weight (number of infections represented by each sequence) based on "proxy_infections"
dat[, 'weight' := proxy_infections / sum(count), by = c('STUSAB', 'yr_wk')]

# Remove NA and INF weights
# index of sequences to be excluded b/c of invalid weights
invalid_weight <- is.na(dat$weight) | is.infinite(dat$weight)

# remove sequences excluded b/c of invalid weights
dat <- subset(x = dat,
              !invalid_weight)


## Survey Design & Weight Trimming ---------------------------------------------

# Define the data date and the week of the data
data_date <- as.Date('2023-05-13')
current_week <- as.numeric(data_date - week0day1) %/% 7

# add the current week to the source data
dat$current_week <- current_week

# create another column for the variants of interest
# this is only used to get (unweighted) counts of the sequences by lineage
dat$VARIANT2 <- as.character(dat$VARIANT)
# group all non-"variants of interest" together
dat[dat$VARIANT %notin% voc, "VARIANT2"] <- "Other"


# specify survey design 
svyDES = survey::svydesign(ids     = ~ SOURCE,
                           strata  = ~ STUSAB + yr_wk,
                           weights = ~ weight,
                           nest    = TRUE, # TRUE = disable checking for duplicate cluster ID's across strata
                           data    = dat)  # not specifying fpc signifies sampling with replacement

# maximum weight; trim weights > 99th percentile
max_weight <- quantile(x = weights(svyDES),
                       probs = .99)

# get the minimum (non-0) weight to replace weights of 0
# weights of 0 can happen if very limited testing data are available for a state
# in a given week. E.g. if 0 of 2 tests were positive, then any sequences that 
# came from that state in that week would represent 0 infections. 
wts <- weights(svyDES)
min_wt <- min(wts[wts > 0])
  
# trim weights using the min & max weights calculated above
svyDES <- survey::trimWeights(design = svyDES,
                              upper = max_weight,
                              lower = min_wt)

# add weights back into the dataframe
dat$weights = weights(svyDES)

## Center "week" for use in Nowcast multinomial model --------------------------

# final date of data to include in model
time_end <- '2023-05-13'
# number of weeks to include in the model
model_weeks <- 21
# maximum week included in the model
model_week_max = as.numeric(as.Date(time_end) - week0day1) %/% 7
# 1 week before first week included in the model
model_week_min = model_week_max - model_weeks
# a midpoint week that will be used to center week values 
model_week_mid = round(model_weeks/2) # center it around 0 instead of using 1:model_weeks
# create a dataframe of old and new values to help visualize things
model_week_df = data.frame(week = 1:model_weeks + model_week_min, # week number (since week0day1)
                           model_week = 1:model_weeks - model_week_mid, # week number for use in nowcast model
                           week_start = (1:model_weeks + model_week_min)*7 + as.Date(week0day1), # date of first day of week
                           week_mid   = (1:model_weeks + model_week_min)*7 + as.Date(week0day1) + 3, # date of midday of week
                           week_end   = (1:model_weeks + model_week_min)*7 + as.Date(week0day1) + 6) # date of final day of week
# add the new model_week info to the dataframe
dat$model_week = dat$week - model_week_min - model_week_mid

# Get the "model_week" of the current week (to be used for making predictions for the current week)
# (this might be the same as the last week in "model_week_df"; it depends on "time_end")
data_week_df = data.frame(
  week       = current_week, 
  model_week = current_week - model_week_min - model_week_mid,
  week_start = (data_date - as.numeric(format(data_date, '%w'))),
  week_mid   = (data_date - as.numeric(format(data_date, '%w'))) + 3,
  week_end   = (data_date - as.numeric(format(data_date, '%w'))) + 6
)

# # a function to convert dates to model_week
# date_to_model_week = function(date){
#   # fractional week (centered on Wednesday)
#   week = as.numeric(as.Date(date) - (week0day1+3)) / 7
#   return(week - model_week_min - model_week_mid)
# }



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Part 3: Calculating Variant Proportions -------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# note: "Proportion" and "Share" are used synonymously in this code

## Weighted variant proportions ------------------------------------------------------

# calculate the variant proportion and confidence interval using survey
# design for:
# - fortnights (HHS regions & nationally)

# get the relevant fortnights from the data
ftnts = sort(unique(dat$FORTNIGHT_END))# create a dataframe with all unique combinations of variants, fortnights, and regions
# create a dataframe with all unique combinations of variants, fortnights, and regions
# (and then reverse column order for convenience)
all.ftnt = expand.grid(Variant = voc,
                       Fortnight_ending = ftnts,
                       USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

# Define a function to calculate the proportion estimates & CI
all_ftnt_ests_function <- function(all.ftnt, svyDES){
  # calculate estimates with 95% CI
  ests <- apply(X = all.ftnt,
                MARGIN = 1,
                FUN = function(rr) myciprop(voc   = rr[3],
                                            geoid = rr[1],
                                            svy = subset(x = svyDES,
                                                          FORTNIGHT_END == rr[2]),
                                            str = FALSE,
                                            level = 0.95))
  return(ests)
} # end all_ftnt_ests_function definition

# calculate weighted proportions (in each fortnight & region) for lineages listed in "voc"
ests <- all_ftnt_ests_function(all.ftnt = all.ftnt, svyDES = svyDES)

# add in the estimates to the dataframe
all.ftnt = cbind(all.ftnt,
                Share    = ests[1,],
                Share_lo = ests[2,],
                Share_hi = ests[3,],
                DF       = ests[4,],
                eff.size = ests[5,],
                cv.mean  = ests[6,],
                deff     = ests[7,])

# make predictions for the "other" variants (following the same steps)
others = expand.grid(Variant          = "Other",
                     Fortnight_ending = ftnts,
                     USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

# function to get the proportion estimates & CI (for "other" variants)
others_ftnt_ests_function <- function(others, svyDES){
  # calculate estimates with 95% CI
  ests.others = apply(X = others,
                      MARGIN = 1,
                      FUN = function(rr) myciprop(voc = voc, # "Other" isn't in the Nowcast model; if voc is a vector, myciprop will return the aggregated proportion. So feed in all the vocs, then subtract from 1 to get "Other" estimates.
                                                  geoid = rr[1],
                                                  svy = subset(svyDES,
                                                              FORTNIGHT_END == rr[2]),
                                                  str = FALSE,
                                                  level = 0.95))

  return(ests.others)
} # end others_ftnt_ests_function definition

# calculate the proportion of "other" variant in each fortnight & region
ests.others <- others_ftnt_ests_function(others = others, svyDES = svyDES)

# add in the estimates to the dataframe (for "other" variants)
others = cbind(others,
               Share    = 1-ests.others[1,],
               Share_lo = 1-ests.others[3,],
               Share_hi = 1-ests.others[2,],
               DF       = ests.others[4,],
               eff.size = ests.others[5,],
               cv.mean  = ests.others[6,],
               deff     = ests.others[7,])

# combine the estimates for the vocs with the estimates for "other" variants
all.ftnt = rbind(all.ftnt,
                 others)

# create a table of counts by variant, time period, and HHS region
raw_counts_REG <- aggregate(count ~ VARIANT2 + FORTNIGHT_END + HHS,
                            data = dat,
                            FUN  = sum,
                            drop = FALSE) # drop = FALSE: keep all combinations, even if no observations

# convert HHS region to character
raw_counts_REG$HHS <- as.character(raw_counts_REG$HHS)

# create a table of counts by variant and time period for the whole US
raw_counts_US <- aggregate(count ~ VARIANT2 + FORTNIGHT_END,
                           data = dat,
                           FUN  = sum,
                           drop = FALSE)

# add a column for HHS region
raw_counts_US <- cbind(raw_counts_US[,1:2],
                       HHS = "USA",
                       count = raw_counts_US[,3])

# combine dataframe of counts by region with dataframe of counts for US
raw_counts <- rbind.data.frame(raw_counts_US,
                               raw_counts_REG)

# merge weighted proportions estimates with sequence counts
all.ftnt2 <- merge(x = all.ftnt,
                   y = raw_counts,
                   by.x = c("USA_or_HHSRegion",
                            "Fortnight_ending",
                            "Variant"),
                   by.y = c("HHS",
                            "FORTNIGHT_END",
                            "VARIANT2"),
                   all.x = T)

# replace NA counts with 0
all.ftnt2[is.na(all.ftnt2$count)==T, "count"] <- 0

#calculate denominator counts by region & time period
dss <- aggregate(count ~ USA_or_HHSRegion + Fortnight_ending,
                 data = all.ftnt2,
                 FUN  = sum)

# change the names of the "count" column
names(dss)[grep("count",names(dss))] <- "denom_count"

# add the denominator counts into the dataframe of results
all.ftnt2 <- merge(x = all.ftnt2,
                   y = dss)

#set the Share 0 and CI limits to NA when the count for a lineage is 0
all.ftnt2$Share = ifelse(test = all.ftnt2$Share != 0 & all.ftnt2$count == 0,
                         yes = 0,
                         no = all.ftnt2$Share)
all.ftnt2$Share_lo = ifelse(test = is.na(all.ftnt2$Share_lo) == F & all.ftnt2$count == 0,
                            yes = NA,
                            no = all.ftnt2$Share_lo)
all.ftnt2$Share_hi = ifelse(test = is.na(all.ftnt2$Share_hi) == F & all.ftnt2$count==0,
                            yes = NA,
                            no = all.ftnt2$Share_hi)

# calculate absolute CI width
all.ftnt2$CI_width = all.ftnt2$Share_hi - all.ftnt2$Share_lo

## generate NCHS flags indicating the reliability of the estimates
# see https://www.cdc.gov/nchs/data/series/sr_02/sr02_175.pdf
# flag estimates with Degrees of Freedom < 8
all.ftnt2$flag_df = as.numeric(all.ftnt2$DF < 8)
# flag estimates with effective size of < 30 (or NA)
all.ftnt2$flag_eff.size = ifelse(test = all.ftnt2$eff.size < 30 |
                                   is.na(all.ftnt2$eff.size) == T,
                                 yes = 1,
                                 no = 0)
# flag estimates with "denominator count" of < 30 (or NA)
# (denominator = count of sequences of all variants in a given region and time period)
all.ftnt2$flag_dss = ifelse(test = all.ftnt2$denom_count < 30 |
                              is.na(all.ftnt2$denom_count) == T,
                            yes = 1,
                            no = 0)
# flag estimates with wide (absolute) confidence intervals
all.ftnt2$flag_abs.ciw = ifelse(test = all.ftnt2$CI_width > 0.30 |
                                  is.na(all.ftnt2$CI_width) == T,
                                yes = 1,
                                no = 0)
# flag estimates with wide (relative) confidence intervals
all.ftnt2$flag_rel.ciw = ifelse(test = ((all.ftnt2$CI_width/all.ftnt2$Share)*100) > 130 |
                                  is.na((all.ftnt2$CI_width/all.ftnt2$Share)*100) == T,
                                yes = 1,
                                no = 0)

# Single identifier for observations that have *any* NCHS flag
all.ftnt2$nchs_flag = ifelse(test = all.ftnt2$flag_df == 1 |
                               all.ftnt2$flag_eff.size == 1 |
                               all.ftnt2$denom_count == 1 |
                               all.ftnt2$flag_abs.ciw == 1 |
                               all.ftnt2$flag_rel.ciw == 1,
                             yes = 1,
                             no = 0)
# Single identifier for observations that have any NCHS flag *other* than the
# degrees of freedom flag.
all.ftnt2$nchs_flag_wodf = ifelse(test = all.ftnt2$flag_eff.size == 1 |
                                    all.ftnt2$denom_count == 1 |
                                    all.ftnt2$flag_abs.ciw == 1 |
                                    all.ftnt2$flag_rel.ciw == 1,
                                  yes = 1,
                                  no = 0)

# write results to file
write.csv(x = all.ftnt2,
          file = 'fortnightly_weighted_proportions.csv',
          row.names = FALSE)


## Nowcast Variant Proportions --------------------------------------------------

### Prep ------------------------------------------------------------------------

# Nowcast model includes all variants listed in "voc" and also the "n_top"
# variants that have the largest weighted variant proportions
n_top = 10 

# define n_top variants over n_recent_weeks
n_recent_weeks = 4

# All variants ordered by (decreasing) weighted variant proportions
us_var = sort(
  prop.table(x = xtabs(formula = weights ~ VARIANT,
                       data = dat, 
                       subset = (dat$week >= current_week - n_recent_weeks))),
  decreasing = TRUE)

# names of all the variants
us_rank = names(us_var)

# variants to include in the nowcast model
# (variants are only included if they were observed during "n_recent_week")
model_vars = us_rank[us_rank %in% c(us_rank[1:n_top], voc)]
# or can just use voc (as long as all voc are in model data)
model_vars = voc[voc %in% dat[ model_week %in% ((1:model_weeks)-model_week_mid)][['VARIANT']] ]

# add variant ranks to dat
# makes sure each seq is assigned a rank/number based on the weighted proportion
# in the last few weeks (1 = most common)
# these are treated as categories and used as the response variable in the 
# multinomial nowcast model
dat$K_US <- match(x = dat$VARIANT, table = model_vars)

# for variants that are not in "model_vars", assign last number, 
# which will correspond to "Other"
dat$K_US[is.na(dat$K_US)] = length(model_vars) + 1

#create a subset of dat that only contains the weeks that will be included
# in multinomial model
moddat = subset(dat,
                model_week %in% ((1:model_weeks)-model_week_mid))

# re-scale model weights 
moddat$wts <- moddat$weights / max(moddat$weights) 
 
#create survey design object based on moddat
mysvy = survey::svydesign(ids     = ~ SOURCE,
                          strata  = ~ STUSAB + yr_wk,
                          weights = ~ wts,
                          nest = TRUE, # TRUE = disable checking for duplicate cluster ID's across strata
                          data = moddat) # not specifying fpc signifies sampling with replacement

### Fit Nowcast models ---------------------------------------------------------
# Fit the nowcast model for regional estimates using the svymultinom function
svymlm_hhs = svymultinom(mod.dat = moddat,
                         mysvy   = mysvy,
                         fmla    = formula("as.numeric(as.factor(K_US)) ~ model_week + as.factor(HHS)"), 
                         model_vars = model_vars)

# Run another model for the entire US
svymlm_us = svymultinom(mod.dat = moddat,
                        mysvy   = mysvy,
                        fmla    = formula("as.numeric(as.factor(K_US))  ~ model_week"), 
                        model_vars = model_vars)

### Aggregate results ----------------------------------------------------------

# generate a matrix that indicates which lineages to aggregate for the nowcast
# Columns are the lineages in the nowcast model
# Rows are the aggregated lineages for output
# a 1 in any cell indicates that the variant indicated by the column should be 
# aggregated into the aggregate indicated by the row.

# Specify the lineages to output (must be a subset of "voc")
output_lineages = model_vars[model_vars %in% voc]

# create an empty aggregation matrix with a column for each voc in the model (plus "OTHER")
agg_var_mat <- matrix(data = 0,
                      nrow = 0,
                      ncol = (length(model_vars)+1))
# specify the column names
colnames(agg_var_mat) <- c(model_vars,"Other")

# get the lineage_expanded for the model_vars
model_vars_expanded <- setNames(voc_lut$lineage_expanded, voc_lut$variant)[ model_vars ]

# potential parent variants are output_lineages_expanded
output_lineages_expanded <- setNames(voc_lut$lineage_expanded, voc_lut$variant)[ output_lineages ]

# get the parent varients for the model_vars_expanded from the output_lineages_expanded
model_var_parents_expanded <- nearest_parent( model_vars_expanded,  output_lineages_expanded )
# get the short names for the parent variants of model_vars
model_var_parents <- setNames( voc_lut$variant, voc_lut$lineage_expanded )[ model_var_parents_expanded ]

# model_var look-up table
model_var_lut <- data.frame(
model_vars = model_vars,
output_model_vars = model_var_parents
)

# only need to convert rows where model_vars != output_model_vars
mvl_sub <- model_var_lut[ model_var_lut$model_vars != model_var_lut$output_model_vars, ]

# get the names of the parent variants that will have sub-lineages aggregated into them
unique_mvl_sub <- unique(mvl_sub$output_model_vars)

# for each unique mvl_sub$output_model_vars, fill in values in the matrix for the variants that will be aggregated into it.
for (i in unique_mvl_sub) {

    # define an extra row to add onto the agg_var_mat
    extra_row <- ifelse( colnames(agg_var_mat) %in%
                            # the variant itself & all the sub-lineages to be aggregated into it
                            c(i, model_var_lut[ model_var_lut$output_model_vars == i, 'model_vars']) ,
                            1,
                            0)

    # add the new row onto the aggregation matrix
    agg_var_mat <- rbind(
        agg_var_mat,
        extra_row
    )
    # give the new row a row name
    row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(i, 'Aggregated')
} # end loop over unique_mvl_sub

# all the variants (not in output_lineages) AND (not aggregated into something else)
# should be aggregated into "Other Aggregated"
other_agg <- base::setdiff(colnames(agg_var_mat)[colSums(agg_var_mat) == 0],
                            output_lineages)
# add the new row (for "Other Aggregated") onto the aggregation matrix
agg_var_mat <- rbind(
    agg_var_mat,
    ifelse(colnames(agg_var_mat) %in% other_agg, 1, 0)
)
# update the row name
row.names(agg_var_mat)[nrow(agg_var_mat)] <- 'Other Aggregated'

# double-check that no variant is aggregated into multiple vocs
if(max(colSums(agg_var_mat)) > 1){
warning(message = paste(
    'Aggregated results are invalid! These variants are being aggregated multiple times:',
    names(agg_var_mat)[colSums(agg_var_mat) > 1],
    '. Fix the aggregation matrix.'))
}

# NOTE! Setting up the agg_var_mat in this way assumes that variants included in
# model_vars but NOT in voc (i.e: setdiff(model_vars, voc) ) will be aggregated into something
# that IS in voc (other than "Other Aggregated"). This is necessary to ensure that the
# "Other Aggregated" row is the SAME for the weighted estimates and the Nowcast estimates
# this is a check to make sure that the assumption stated above is being met.
sub_mat <- agg_var_mat[row.names(agg_var_mat) != 'Other Aggregated',setdiff(model_vars, voc), drop = FALSE]
if(!all(colSums(sub_mat) > 0)){
    problem_vocs <- colnames(sub_mat)[colSums(sub_mat) == 0]
    warning(message = paste0('Not all variants in "model_vars" are aggregated into something in "voc". ',
                                paste(problem_vocs, collapse = ', '),
                                ' is in model_vars but is not aggregated into anything in voc. Therefore Nowcast & weighted estimates will differ!'))
}

# print a warning if any columns have totals > 1
if(any(colSums(agg_var_mat)>1)) warning(paste0('agg_var_mat not correctly specified. Some variants are aggregated more than once.', agg_var_mat))


### Fortnightly estimates -----
#define fortnights and regions to get nowcasts for
# use the final 2 weeks of observed data
proj_ftnts = as.Date(tail(ftnts, 2))
# and add on 2 fortnights into the future
proj_ftnts = sort(unique(c(proj_ftnts,
                            proj_ftnts + 14,
                            proj_ftnts + 28)))

# create a dataframe with all regions for predictions
dfs = expand.grid(USA_or_HHSRegion = c("USA", as.character(1:10)))

# create an empty object to hold predicted values for each region
proj.res = c()

# cycle over regions
for (rgn in dfs$USA_or_HHSRegion){
  # cycle over time period
  for (ftn in proj_ftnts) {

    # get the model fit object & geiod
    if (rgn=="USA") {
      mlm   = svymlm_us$mlm
      geoid = rgn
    } else {
      mlm   = svymlm_hhs$mlm
      geoid = as.numeric(rgn)
    }

    # get the timepoint to use for predictions for each fortnight
    # (use the midpoint of each fortnight, and convert to the same unit of time as is used in the model)
    wk = as.numeric(as.Date(ftn, origin="1970-01-01") - (week0day1+3) - 6.5) / 7
    # convert week to model_week
    wk = wk - model_week_min - model_week_mid
    
    # get the estimates (and SE) for the given place & time
    ests = se.multinom(mlm = mlm,
                        newdata_1row = data.frame(
                          model_week = wk,
                          HHS = geoid
                        ),
                        composite_variant = agg_var_mat)

    # calculate the growth rate
    #   "growth rate" here is the derivative (i.e. slope) of the predicted proportion
    #   hence, it is an instantaneous rate of change. 
    #   multinomial model estimated proportion of variant i at time t:  p_i(t) = exp(b_0i + b_1i * t) / sum_j(exp(b_0j + b_1j * t)), where j is all the variants in the model
    #   the derivative of the log of this function is: dlog(p_i(t))/dt = b_1i - sum_j(p_j(t) * b_1j)
    #   then we convert back to the response scale by taking the exponent, 
    #   convert from a fraction to a percent by multiplying by 100
    #   and subtract 100 so that "no change" is 0 instead of (100 percent of current value)
    #   A growth rate of 50 means that the instantaneous rate of change of the 
    #   estimated proportion is a 50% increase over 1 week (the unit of time in the model is week).
    #   The growth rates are not constant through time. A variant with a proportion of 50% and
    #   a growth rate of 50% will not reach 75% after one week. It would if the instantaneous 
    #   growth rate were maintained for the full week, but growth rates typically decline with time.
    gr = with(ests,
              100 * exp(b_i - sum(p_i * b_i)) - 100)

    # calculate the SE of the growth rate
    se.gr = with(data = ests,
                  expr = 100 * exp(sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)
    # p_i    = predicted probability (proportional representation/frequency)
    # se.p_i = SE of predicted frequency
    # b_i    = coefficient value (for intercept & week)
    # se.b_i = SE of coefficient value

    # calculate growth rates without rescaling (to use for doubling times & confidence intervals)
    se.gr_link = with(data = ests,
                      expr = sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))
    gr_link = with(data = ests,
                  expr = (b_i - sum(p_i * b_i)))
    gr_lo_link = gr_link - 1.96 * se.gr_link
    gr_hi_link = gr_link + 1.96 * se.gr_link

    gr_lo = 100 * exp(gr_lo_link) - 100
    gr_hi = 100 * exp(gr_hi_link) - 100

    # calculate doubling time (multiply by 7 to convert weeks to days)
    doubling_time    = log(2)/gr_link * 7
    doubling_time_lo = log(2)/gr_lo_link * 7
    doubling_time_hi = log(2)/gr_hi_link * 7

    # empty dataframe to store growth rates (and doubling times) for aggregated variants
    gr_agg <- data.frame( variant = rownames(agg_var_mat),
                          gr    = NA,
                          se.gr = NA,
                          gr_lo = NA,
                          gr_hi = NA,
                          dt    = NA,
                          dt_lo = NA,
                          dt_hi = NA)
    # extract the growth rates for the aggregated variants
    for(r in 1:nrow(agg_var_mat)){
      # if nothing is actually being aggregated, just get the growth rate of the individual component
      if(unname(rowSums(agg_var_mat)[r]) == 1){
        col_ind <- which(agg_var_mat[r,]>0)
        gr_agg[r,'gr']    <- gr[col_ind]
        gr_agg[r,'se.gr'] <- se.gr[col_ind]
        gr_agg[r,'gr_lo'] <- gr_lo[col_ind]
        gr_agg[r,'gr_hi'] <- gr_hi[col_ind]
        gr_agg[r,'dt']    <- doubling_time[col_ind]
        gr_agg[r,'dt_lo'] <- doubling_time_lo[col_ind]
        gr_agg[r,'dt_hi'] <- doubling_time_hi[col_ind]
      } else {
        # if there are component variants, then take the weighted mean (weighted by estimated proportion) to get the aggregated growth rate
        col_ind <- unname(which(agg_var_mat[r,]>0))
        gr_agg[r,'gr']    <- sum(   gr[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
        gr_agg[r,'gr_lo'] <- sum(gr_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
        gr_agg[r,'gr_hi'] <- sum(gr_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
        gr_agg[r,'dt']    <- sum(   doubling_time[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
        gr_agg[r,'dt_lo'] <- sum(doubling_time_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
        gr_agg[r,'dt_hi'] <- sum(doubling_time_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
      }
    }

    # format estimates into dataframe with relevant info
    ests = data.table::data.table(
      # region
      USA_or_HHSRegion = rgn,
      # fortnight
      Fortnight_ending = as.Date(ftn, origin="1970-01-01"),
      # all the variants
      Variant = c(model_vars,
                  "Other",
                  row.names(ests$composite_variant$matrix)),
      # all the estimated proportions
      Share = c(ests$p_i,
                ests$composite_variant$p_i),
      # SE of estimated proportions
      se.Share = c(ests$se.p_i,
                  ests$composite_variant$se.p_i),
      # all the estimated growth rates (and CI)
      growth_rate    = c(gr,    gr_agg$gr),
      growth_rate_lo = c(gr_lo, gr_agg$gr_lo),
      growth_rate_hi = c(gr_hi, gr_agg$gr_hi),
      # all the estimated doubling times (and CI)
      doubling_time    = c(doubling_time,    gr_agg$dt),
      doubling_time_lo = c(doubling_time_lo, gr_agg$dt_lo),
      doubling_time_hi = c(doubling_time_hi, gr_agg$dt_hi)
    )

    # Get binomial CI using the estimated proportion & SE
    binom.ci = apply(X = ests,
                    MARGIN = 1,
                    FUN = function(rr) svyCI(p = as.numeric(rr[4]),
                                             s = as.numeric(rr[5]), 
                                             conf.level = 0.95))

    # add the CI into the estimates dataframe
    ests$Share_lo = binom.ci[1,]
    ests$Share_hi = binom.ci[2,]

    # add the estimates for this specific place & time to the results
    proj.res = rbind(proj.res,
                    ests)
  } # end loop over fortnights
} # end loop over regions

# Extract and save the results for the aggregated variants (as specified in agg_var_mat)
# exclude variants that have been aggregated into other groups (i.e. have value > 0 in the matrix)
agg_lineages <- colnames(agg_var_mat)[colSums(agg_var_mat)>0]
if("Other" %notin% agg_lineages) agg_lineages <- c(agg_lineages, "Other Aggregated")
# get the results that were not aggregated into something else
results_agg = proj.res[Variant %notin% agg_lineages]

# double-check that the proportions add up to 1 each time period
if (!all(results_agg[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'Fortnight_ending')][,unique(round(total_share, 5))] == 1)){
  warning('The total proportion does not add up to 100% in each time period!')
} else {
  # save the results to file
  write.csv(x = results_agg,
            file = 'fortnightly_nowcast_proportions_aggregated.csv',
            row.names = FALSE)
}

# Extract and save the results for all model_vars (i.e. exclude the aggregated variants)
# exclude the lineages that were aggregated (other than "Other")
drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]

# Only include variants that are NOT in the list provided
results_nonagg = proj.res[Variant %notin% c(drop_lin, "Other", colnames(agg_var_mat['Other Aggregated',,drop=F])[agg_var_mat['Other Aggregated',,drop=F] > 0])]

# Double-check that the shares add up to 1 each time period
if (!all(results_nonagg[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'Fortnight_ending')][,unique(round(total_share, 5))] == 1)){
  warning('The total proportion does not add up to 100% in each time period!')
} else {
  write.csv(x = results_nonagg,
            file = 'fortnightly_nowcast_proportions.csv',
            row.names = FALSE)
}
