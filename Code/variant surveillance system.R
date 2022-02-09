# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Code accompanying MMWR on Genomic Surveillance 
# Date: 2021-05-11; updated 2022-01-18
# Author: Prabasaj Paul
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Part 1: Data Prep & weight calculation ---------------------------------------
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

# define "%notin%" function
`%notin%` <- Negate(`%in%`)

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
dat = data.table::data.table(dat)



## Some general parameters -----------------------------------------------------

# starting date for defining weeks
week0day1 = as.Date("2020-01-05")

# Variants to be included in the analysis (not limited to VOCs; can be any Pango lineage)
voc = c("B.1.1.7",   # Alpha
        "P.1",       # Gamma
        "B.1.617.2", # Delta
        "C.37",      # Lambda
        "B.1.621",   # Mu 
        "B.1.1.529") # Omicron

# fewer variants runs faster
voc = c('B.1.617.2', 
        'B.1.1.529', 
        'AY.103')

## Aggregate sublineages -------------------------------------------------------
# all the AY sublineages
AY = sort(unique(dat$VARIANT)[grep("AY",unique(dat$VARIANT))])
# just the AY sublineages that are to be aggregated (i.e. not listed in "voc")
AY = AY[which(AY %notin% voc)]
# aggregate sublineages to the parent lineage
dat[dat$VARIANT %in% AY, "VARIANT"] <- "B.1.617.2"

P1 = sort(unique(dat$VARIANT)[grep("P\\.1.",unique(dat$VARIANT))])
P1 = P1[which(P1 %notin% voc)] #vector of the P1s to aggregate
dat[dat$VARIANT %in% P1, "VARIANT"] <- "P.1"

Q = sort(unique(dat$VARIANT)[grep("Q\\.",unique(dat$VARIANT))])
Q = Q[which(Q %notin% voc)] #vector of the Qs to aggregate
dat[dat$VARIANT %in% Q, "VARIANT"] <- "B.1.1.7"

B351 = sort(unique(dat$VARIANT)[grep("B\\.1\\.351\\.",unique(dat$VARIANT))])
B351 = B351[which(B351 %notin% voc)] #vector of the B351s to aggregate
dat[dat$VARIANT %in% B351,"VARIANT"] <- "B.1.351"

B621 = sort(unique(dat$VARIANT)[grep("B\\.1\\.621\\.",unique(dat$VARIANT))])
B621 = B621[which(B621 %notin% voc)] #vector of the B621s to aggregate
dat[dat$VARIANT %in% B621,"VARIANT"] <- 'B.1.621'

B429 = sort(unique(dat$VARIANT)[grep("B\\.1\\.429",unique(dat$VARIANT))])
B429 = B429[which(B429 %notin% voc)] #vector of the B429s to aggregate
dat[dat$VARIANT %in% B429,"VARIANT"] <- "B.1.427"

B529 = sort(unique(dat$VARIANT)[grep("(B\\.1\\.1\\.529)|(BA\\.[0-9])",unique(dat$VARIANT))])
B529 = B529[which(B529 %notin% voc)] #vector of the B529s to aggregate
dat[dat$VARIANT %in% B529,"VARIANT"] <- "B.1.1.529"

## SGTF over-sampling weights --------------------------------------------------
## weighting adjusted for oversampling/targeted sampling of SGTF+ specimens sent for sequencing
# counts of lineages by contractor name & SGTF upsampling status
sgtf.1 = table(dat$VARIANT,
               paste(dat$SOURCE,
                     dat$SGTF_UPSAMPLING))

# Get the labs/columns with s-gene target failure oversampling
# (identified by having more SGTF_UPSAMPLING==TRUE tests than SOURCE=="17" tests)
sgtf.vars = rownames(sgtf.1)[
  sgtf.1[, grep(pattern = "TRUE$", x = colnames(sgtf.1))] > 
    sgtf.1[, grep(pattern = "17 FALSE$", x = colnames(sgtf.1))] ]

# calculate SGTF weights using sets of states and weeks where targeted samples were sequenced
sgtf.sub = unique(
  subset(
    x = dat,
    SGTF_UPSAMPLING
  )[, c("STUSAB", "yr_wk")]
)

# fit a logistic model to estimate probability of a sequence having been 
# sequenced b/c of SGTF
sgtf.glm = glm(
  I(VARIANT %in% sgtf.vars) ~ I(SOURCE %in% "17") + STUSAB + yr_wk,
  family="binomial",
  subset(dat, yr_wk %in% sgtf.sub$yr_wk & 
           STUSAB %in% sgtf.sub$STUSAB))

# add SGTF upsampling weights to the subset data
sgtf.sub$sgtf_weights = predict(
  object = sgtf.glm, 
  newdata = cbind(SOURCE = c("17"), 
                  sgtf.sub), 
  type="response"
) / predict( 
  object = sgtf.glm, 
  newdata = cbind(SOURCE = c("OTHER"), 
                  sgtf.sub), 
  type="response"
)

# merge SGTF upsampling weights into the main dataset
dat <- merge(x = dat, 
             y = sgtf.sub, 
             by = c('STUSAB', 'yr_wk'), 
             all.x = TRUE)

# replace NA weights with 1
dat[ is.na(dat$sgtf_weights), 'sgtf_weights'] <- 1


## Survey Weights --------------------------------------------------------------
# Estimation of the infection weights involves estimating the (unobserved) number
# of infections from test results. There is no reliable and precise method for this
# yet, so these weights are subject to considerable uncertainties.
# We use methods in:
#     https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm
#     https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
# Other strategies are also available; e.g.
#     https://covid19-projections.com/estimating-true-infections-revisited/


# calculate simple adjusted weights
dat[, 
    "SAW" := sqrt(state_population/TOTAL) * POSITIVE/sum(1/sgtf_weights)/sgtf_weights,
    .(STUSAB, yr_wk)]
# calculate simple adjusted weights for states that lack testing data
# (use average values from HHS region)
dat[, 
    "SAW_ALT" := state_population*HHS_INCIDENCE / sum(1/sgtf_weights)/sgtf_weights,
    .(STUSAB, yr_wk)]
# use "SAW_ALT" when SAW is NA
dat[, "SIMPLE_ADJ_WT" := ifelse(is.na(SAW), SAW_ALT, SAW)]
# remove "SAW" and "SAW_ALT" columns
dat[, c("SAW", "SAW_ALT") := .(NULL, NULL)]


# Make sure all survey weights are valid numbers
# (this shouldn't be necessary, just precautionary)
dat = subset(x = dat,
             !is.na(SIMPLE_ADJ_WT) &
               SIMPLE_ADJ_WT < Inf)

## Survey Design & Weight Trimming ---------------------------------------------

# weights are trimmed to avoid overly influential samples
# current week 
current_week = as.numeric(Sys.Date() - week0day1) %/% 7

# add the current week to the source data
dat$current_week = current_week

# create another column for the variants of interest
# this is only used to get (unweighted) counts of the sequences by lineage
dat$VARIANT2 = as.character(dat$VARIANT)
# group all non-"variants of interest" together
dat[dat$VARIANT %notin% voc, "VARIANT2"] <- "Other"


# specify survey design 
svyDES = survey::svydesign(ids     = ~ SOURCE,
                           strata  = ~ STUSAB + yr_wk,
                           weights = ~ SIMPLE_ADJ_WT,
                           nest    = TRUE, # TRUE = disable checking for duplicate cluster ID's across strata
                           data    = dat)  # not specifying fpc signifies sampling with replacement

# maximum weight; trim weights > 99th percentile
max_weight <- quantile(x = weights(svyDES),
                       probs = .99)

# get the minimum weight to replace weights of 0
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
time_end <- '2022-01-01'
# number of weeks to include in the model
model_weeks <- 20
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
                        ...)
    )
  } else {
    # use a modified version of svyciprop to limit effective sample size so that
    # it's never > observed sampled size.
    res = suppressWarnings(
      svycipropkg(~VOC,
                  design = srv_design,
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

# function to calculate & format a date based on # of weeks from datadate
week_label = function(week_from_current,
                      datadate = Sys.Date()) {
  # returns the date of the first day of a future week in "%m-%d" format
  
  # the current day of the week
  day_of_week = as.numeric(strftime(as.Date(datadate), format="%w"))
  
  # start of current week
  start_of_week = as.Date(datadate) - day_of_week
  
  # formatted date of a future week
  strftime(x = start_of_week + 7 * week_from_current,
           format = "%m-%d")
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
                       fmla = formula("as.numeric(as.factor(K_US)) ~ week + as.factor(HHS)"),
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
  
  # add the variant names into the results to make the results easier to read
  # (first variant is used as comparison for all others)
  # names(rval$estimates) <- modvars$Variant[-1]
  
  # use a "try-catch" function to try to calculate the inverse of the negative  
  # Hessian, which is used below.  
  invinf <- tryCatch( 
    { 
      solve(-multinom_geoid$Hessian) 
    }, 
    error = function(cond) { 
      return(NA) 
    } 
  ) 

  # if the Hessian is not invertible, avoid errors by not running the rest of the function
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
    # if the Hessian is solvable, proceed:
    
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
                       week,
                       geoid = "USA",
                       composite_variant = NULL) {
  # Arguments:
  #  ~  mlm:   model output, with Hessian;
  #  ~  week:  focal week (model_week) to calculate the SE for
  #  ~  geoid: geographic region; options include "USA" or HHS Region (1:10)
  #  ~  composite_variant: (if not NA) is a matrix with one column per model
  #                        lineage and one row per composite variant. Matrix
  #                        element of 1 marks each component lineage (column)
  #                        for each composite variant (row).
  #                        Example: matrix(c(1, 0, 0, 1, 0, 0), nrow = 1)
  #                        designates a variant comprising the first and fourth
  #                        lineages in the model estimates for composites will
  #                        be appended at end of p_i and se.p_i
  
  # If mlm is without Hessian, all se's are set to zero
  # mlm not geographically stratified for geoid="USA": ~ (week - current_week)
  # mlm, for HHS Regions: ~ week + HHS (1 as reference level)
  
  # Output: list of 5 items:
  #  ~  p_i:     predicted values for each clade/variant (numeric vector)
  #  ~  se.p_i:  SE of predicted values (numeric vector)
  #  ~  b_i:     coefficient (beta) values (numeric vector)
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
      vc = rep(0, length(cf)^2)
    }
  }
  
  # Rearrange terms to ease pulling out relevant submatrix
  # convert symmetrical matrix with rows & columns for each coefficient:
  # 4-d array with dim 1 for each covariate
  #                dim 2 for each outcome class
  #                dim 3 for each covariate
  #                dim 4 for each outcome class
  dim(vc) = rep(x = rev(dim(cf)),
                times = 2)
  
  # Boolean: is the geoid the reference level
  ref_geoid = ((length(mlm$xlevels) == 0) || (geoid == mlm$xlevels[[1]][1]))
  
  # get the number/index of the geoid/area being modeled
  if (ref_geoid) { # if the geoid is the reference level
    n_geoid = NULL
  } else { # if it's not the reference level
    n_geoid = which(mlm$xlevels[[1]] == geoid)
  }
  
  # Set value for "hhs1tail"
  # index of column identifying column for the hhs region in question
  if (ref_geoid) { # if the geoid is the reference level
    hhs1tail = NULL
  } else { # if it's not the reference level
    hhs1tail = n_geoid + 1
  }
  
  # save values to vector
  # only want the coefficients for the intercept, the week in question, and the
  # hhs region in question
  indices = c(1, 2, hhs1tail)
  
  # reset value for "hhs1tail"
  # indicator (0/1) covariate value (1 for the hhs region in question)
  if (ref_geoid) {
    hhs1tail = NULL
  } else {
    hhs1tail = 1
  }
  
  # save values to another vector
  coeffs = c(1, week, hhs1tail)
  
  # b is vector of coefficients of time (week)
  # subset the variance-covariance matrix to the relevant parameters
  sub.vc = vc[indices,,indices,]
  
  # get the linear predictor values
  y_i = c(0, cf[, indices] %*% coeffs)
  
  # get the variance-covariance matrix of the linear predictor
  # (after adjusting for survey design)
  vc.y_i = outer(1:dim(sub.vc)[2],
                 1:dim(sub.vc)[4],
                 Vectorize(function(i, j) c(coeffs %*% sub.vc[, i, , j] %*% coeffs)))
  
  # add in row and column with 0
  vc.y_i = rbind(0, cbind(0, vc.y_i))
  
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
  
  # output
  list(p_i = p_i, # predicted values
       se.p_i = se.p_i, # SE of predicted values
       b_i = c(0, cf[, 2]), # coefficient value for week
       se.b_i = c(0, sqrt(diag(vc[2,,2,]))), # SE of the coefficient
       composite_variant = composite_variant) # predicted probabilities (and SE) for the composite variants
}

# Function to get binomial confidence interval based on the point estimated
# proportion and associated SE
# (based on output from svymultinom & se.multinom functions)
svyCI = function(p, s) {
  # p = point estimate
  # s = estimated standard error
  
  # if se is 0, n will be Inf (possibly because of non-invertible Hessian); return CI of 0,0
  if(s == 0){
    return(c(0,0))
  } else {
    # calculate the sample size
    n = p*(1-p)/s^2
    
    # calculate the confidence interval using a proportion test
    out = prop.test(x = n * p,   # a vector of counts of successes
                    n = n        # a vector of counts of trials
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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Part 2: Calculating Variant Proportions --------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# note: "Proportion" and "Share" are used synonymously in this code

## Weighted variant proportions ------------------------------------------------------

# calculate the variant proportion and confidence interval using survey
# design for:
# - weeks  (HHS regions & nationally)

# subset data to only include data since 2 May, 2021
# and (end of week) older than "time_end"
dat2 <- dat[yr_wk >= "2021-05-02" & yr_wk <= (as.Date(time_end) - 6), ]
    
# add a column for the date of the final day of each week
dat2$WEEK_END = as.Date(dat2$yr_wk) + 6
    
# get the unique weeks that are in the appropriate time frame
wks = sort(unique(dat2$yr_wk))

# create a dataframe with variant, time period, and region
# (and reorder columns for convenience)
all.wkly = expand.grid(Variant = voc,
                       Week_of = wks,
                       USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

# get predictions & CI for each variant in each week and region
ests = apply(X = all.wkly,
             MARGIN = 1,
             FUN = function(rr) myciprop(voc   = rr[3],
                                         geoid = rr[1],
                                         svy   = subset(svyDES, yr_wk == rr[2]),
                                         str   = FALSE))

# add the predictions into the dataframe
all.wkly = cbind(all.wkly,
                 Share    = ests[1,],
                 Share_lo = ests[2,],
                 Share_hi = ests[3,],
                 DF       = ests[4,],
                 eff.size = ests[5,],
                 cv.mean  = ests[6,])

# get predictions for the non-focal ("other") variants grouped together
# (and reorder columns for convenience)
others = expand.grid(Variant = "Other",
                     Week_of = wks,
                     USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

# get predictions & CI for "other" (grouped) variant in each week and region
ests.others = apply(X = others,
                    MARGIN = 1,
                    FUN = function(rr) myciprop(voc = voc,
                                                geoid = rr[1],
                                                svy = subset(svyDES, yr_wk == rr[2]),
                                                str =  FALSE))

# add the predictions into the dataframe
others = cbind(others,
               Share    = 1-ests.others[1,],
               Share_lo = 1-ests.others[3,],
               Share_hi = 1-ests.others[2,],
               DF       = ests.others[4,],
               eff.size = ests.others[5,],
               cv.mean  = ests.others[6,])

# combine estimates for individual variants with the estimates for the "other"
# (grouped) variants
all.wkly = rbind(all.wkly,
                 others)

# Add a column for the last day of each week
all.wkly$WEEK_END = as.Date(all.wkly$Week_of) + 6


## generate sequence counts by lineage, location and date
# raw counts of each variant in each week and each region
raw_counts_REG <- aggregate(count ~ VARIANT2 + WEEK_END + HHS,
                            data = dat2,
                            FUN  = sum,
                            drop = FALSE)

# convert region names to character
raw_counts_REG$HHS <- as.character(raw_counts_REG$HHS)

# raw counts for each variant in each week nationally
raw_counts_US <- aggregate(count ~ VARIANT2 + WEEK_END,
                           data = dat2,
                           FUN  = sum,
                           drop = FALSE)

# add in a column for HHS region = "USA"
raw_counts_US <- cbind(raw_counts_US[,1:2],
                       HHS   = "USA",
                       count = raw_counts_US$count)

# combine the regional counts with the national counts
raw_counts <- rbind.data.frame(raw_counts_US,
                               raw_counts_REG)

# merge sequence counts with weighted proportions estimates
all.wkly2 <- merge(x = all.wkly,
                   y = raw_counts,
                   by.x = c("USA_or_HHSRegion",
                            "WEEK_END",
                            "Variant"),
                   by.y = c("HHS",
                            "WEEK_END",
                            "VARIANT2"),
                   all = T)

# Set counts to 0 if current count is NA and variant != "Other"
all.wkly2[is.na(all.wkly2$count)==T ,"count"] <- 0

# calculate denominator counts
# (i.e. raw counts per region & week)
dss <- aggregate(count ~ USA_or_HHSRegion + WEEK_END,
                 data = all.wkly2,
                 FUN  = sum)

# change name of denominator counts
names(dss)[grep("count",names(dss))] <- "denom_count"

# add denominator counts into the dataframe of raw counts
all.wkly2 <- merge(x = all.wkly2,
                   y = dss)

#calculate absolute CI width
all.wkly2$CI_width = all.wkly2$Share_hi - all.wkly2$Share_lo

#set the Share (i.e. Proportion) to 0 and CI limits to NA when the count for a lineage is 0
all.wkly2$Share=ifelse(all.wkly2$Share!=0 & all.wkly2$count==0,0,all.wkly2$Share)
all.wkly2$Share_lo=ifelse(is.na(all.wkly2$Share_lo)==F & all.wkly2$count==0,NA,all.wkly2$Share_lo)
all.wkly2$Share_hi=ifelse(is.na(all.wkly2$Share_hi)==F & all.wkly2$count==0,NA,all.wkly2$Share_hi)

## generate NCHS flags indicating the reliability of the estimates
# see https://www.cdc.gov/nchs/data/series/sr_02/sr02_175.pdf
# flag estimates with Degrees of Freedom < 8
all.wkly2$flag_df = ifelse(test = all.wkly2$DF < 8,
                           yes = 1,
                           no = 0)
# flag estimates with effective size of < 30 (or NA)
all.wkly2$flag_eff.size = ifelse(test = all.wkly2$eff.size < 30 |
                                   is.na(all.wkly2$eff.size) == T,
                                 yes = 1,
                                 no = 0)
# flag estimates with "denominator count" of < 30 (or NA)
# (denominator = count of sequences of all variants in a given region and time period)
all.wkly2$flag_dss = ifelse(test = all.wkly2$denom_count < 30|
                              is.na(all.wkly2$denom_count) == T,
                            yes = 1,
                            no = 0)
# flag estimates with wide (absolute) confidence intervals
all.wkly2$flag_abs.ciw = ifelse(test = all.wkly2$CI_width > 0.30 |
                                  is.na(all.wkly2$CI_width) == T,
                                yes = 1,
                                no = 0)
# flag estimates with wide (relative) confidence intervals
all.wkly2$flag_rel.ciw = ifelse(test = ((all.wkly2$CI_width/all.wkly2$Share)*100) > 130 |
                                  is.na((all.wkly2$CI_width/all.wkly2$Share)*100) == T,
                                yes = 1,
                                no = 0)

# Single identifier for observations that have *any* NCHS flag
all.wkly2$nchs_flag = ifelse(test = all.wkly2$flag_df == 1 |
                               all.wkly2$flag_eff.size == 1 |
                               all.wkly2$denom_count == 1 |
                               all.wkly2$flag_abs.ciw == 1 |
                               all.wkly2$flag_rel.ciw == 1,
                             yes = 1,
                             no = 0)
# Single identifier for observations that have any NCHS flag *other* than the
# degrees of freedom flag.
all.wkly2$nchs_flag_wodf = ifelse(all.wkly2$flag_eff.size == 1 |
                                    all.wkly2$denom_count == 1 |
                                    all.wkly2$flag_abs.ciw == 1 |
                                    all.wkly2$flag_rel.ciw == 1,
                                  yes = 1,
                                  no = 0)

# select the columns to save
all.wkly2 = all.wkly2[,c("USA_or_HHSRegion",
                         "WEEK_END",
                         "Variant",
                         "Share",
                         "Share_lo",
                         "Share_hi",
                         "count",
                         "denom_count",
                         "DF",
                         "eff.size",
                         "CI_width",
                         "nchs_flag",
                         "nchs_flag_wodf")]

# sort by HHS region
all.wkly2 <- all.wkly2[order(all.wkly2$USA_or_HHSRegion),]

# save the results to file
write.csv(
  x = all.wkly2,
  file = "variant_share_weekly_weighted.csv",
  row.names = FALSE
)

## Nowcast Variant Proportions --------------------------------------------------

### Prep ------------------------------------------------------------------------

# Model is fit to n_top variants that have the largest weighted variant proportions
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
# Can use those in n_top or in voc
model_vars = us_rank[us_rank %in% c(us_rank[1:n_top], voc)]
# or can just use voc
model_vars = voc

# add variant ranks to dat
# makes sure each seq is assigned a rank/number based on the weighted proportion
# in the last few weeks (1 = most common)
# these are treated as categories and used as the response variable in the 
# multinomial nowcast model
dat$K_US <- match(x = dat$VARIANT, table = model_vars)

# for variants that are not in "model_vars", assign last number, 
# which will correspond to "Other"
dat$K_US[is.na(dat$K_US)] = length(model_vars) + 1

#create a subset of src.dat that only contains the weeks that will be included
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
# Fit the nowcast model using the svymultinom function
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

# optionally aggregate Delta's AY sublineages for results
# Check to see which AY lineages are in model_vars
AY_agg = model_vars[grep("AY", model_vars)]

# all variants to be aggregated into the "other" category
Other_agg = model_vars[model_vars %notin% voc]

# generate a matrix that indicates which lineages to aggregate for the nowcast
# Columns are the lineages in the nowcast model, so all the defined lineages
#  plus the "other" lineage
# Rows are the aggregated lineages wanted
agg_var_mat <- matrix(data = 0,
                      nrow = 2,
                      ncol = (length(model_vars)+1))
colnames(agg_var_mat) <- c(model_vars,"Other")

# Fill in matrix values: if lineage is to be aggregated to parent lineage in
# given row, then value = 1, else value = 0
agg_var_mat[1,] <- ifelse(colnames(agg_var_mat) %in% c("B.1.617.2", AY_agg),1,0)
agg_var_mat[2,] <- ifelse(colnames(agg_var_mat) %in% c(Other_agg, "Other"),1,0)
row.names(agg_var_mat) <-c("Delta Aggregated", "Other Aggregated")

# dates for which to make predictions using the nowcast model
prediction_dates = seq.Date(from = as.Date('2021-11-01'),
                            to = as.Date('2022-01-15'),
                            by = 1)

# create an empty object to hold predicted values for each region
proj.res <- c()

# cycle over regions
for (rgn in c('USA', 1:10)){
  # cycle over dates
  for (pd in prediction_dates) {
    
    # get the model fit & geoid
    if (rgn=="USA") {
      mlm   = svymlm_us$mlm
      geoid = rgn
    } else {
      mlm   = svymlm_hhs$mlm
      geoid = as.numeric(rgn)
    }
    
    # get the week for the given timepoint
    wk_date = as.Date(pd, origin="1970-01-01")
    wk_num  = as.numeric(wk_date - week0day1) / 7
    # convert to model_week
    wk      = wk_num - model_week_min - model_week_mid
    
    # get the estimates (and SE) for the given place & time
    ests = se.multinom(mlm = mlm,
                       week = wk,
                       geoid = geoid,
                       composite_variant = agg_var_mat)
    
    # calculate the SE of the estimated growth rate
    se.gr = with(data = ests,
                 expr = 100 * exp(sqrt(se.b_i^2 + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)
    # p_i    = predicted probability (proportional representation/frequency)
    # se.p_i = SE of predicted frequency
    # b_i    = coefficient value (for intercept & week)
    # se.b_i = SE of coefficient value
    
    # calculate the estimated growth rate
    # here "growth rate" = week-over-week growth of log odds
    gr = with(data = ests,
              expr = 100 * exp(b_i - sum(p_i * b_i)) - 100)
    
    # Get approximate 95% confidence intervals for growth rate
    se.gr_link = with(data = ests,
                      expr = sqrt(se.b_i^2 + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))
    gr_link = with(data = ests,
                   expr = (b_i - sum(p_i * b_i)))
    gr_lo_link = gr_link - 1.96 * se.gr_link
    gr_hi_link = gr_link + 1.96 * se.gr_link
    gr_lo = 100 * exp(gr_lo_link) - 100
    gr_hi = 100 * exp(gr_hi_link) - 100
      
    # calculate doubling time (doubling time of log odds in days)
    doubling_time    = log(2)/gr_link * 7
    doubling_time_lo = log(2)/gr_lo_link * 7
    doubling_time_hi = log(2)/gr_hi_link * 7
    
    
    # format estimates into dataframe with relevant info
    ests_df = data.frame(
      USA_or_HHSRegion = rgn,
      DATE = as.Date(pd,
                     origin="1970-01-01"),
      Variant = c(model_vars,
                  "Other",
                  row.names(ests$composite_variant$matrix)),
      Share = c(ests$p_i,
                ests$composite_variant$p_i),
      se.Share = c(ests$se.p_i,
                   ests$composite_variant$se.p_i),
      growth_rate = c(gr,
                      rep(NA, length(ests$composite_variant$p_i))),
      se.growth_rate = c(se.gr,
                         rep(NA, length(ests$composite_variant$p_i))),
      growth_rate_lo = c(gr_lo,
                         rep(NA, length(ests$composite_variant$p_i))),
      growth_rate_hi = c(gr_hi,
                         rep(NA, length(ests$composite_variant$p_i))),
      doubling_time = c(doubling_time,
                        rep(NA, length(ests$composite_variant$p_i))),
      doubling_time_lo = c(doubling_time_lo,
                           rep(NA, length(ests$composite_variant$p_i))),
      doubling_time_hi = c(doubling_time_hi,
                           rep(NA, length(ests$composite_variant$p_i)))
    )
    
    # Get binomial CI from p_i and se.p_i
    binom.ci = apply(X = ests_df,
                     MARGIN = 1 ,
                     FUN = function(rr) svyCI(p = as.numeric(rr[4]),
                                              s = as.numeric(rr[5])))
    # add the CI into the estimates dataframe
    ests_df$Share_lo = binom.ci[1,]
    ests_df$Share_hi = binom.ci[2,]
    
    # add the estimates for this specific place & time to the results
    proj.res = rbind(proj.res,
                     ests_df)
  } # end loop over weeks
} # end loop over regions

# In example code, "Delta Aggregated" = B.1.617.2 and AY.103 combined
# save results
write.csv(x = proj.res,
          file = "nowcast_weekly.csv",
          row.names = FALSE)


## Jurisdictional estimates ----------------------------------------------------
# this is very similar to the Weighted Variant Share estimates above, but instead
# of producing estimates for each week, estimates are made for rolling 4-week bins
# b/c sample sizes from individual state-week combinations are often too small
# to produce estimates for every week for every state

# Last date of data that gets included in jurisdictional estimates
# typically exclude the most recent 2 weeks of data b/c of delay in sequencing & 
# reporting samples
current_day <- as.Date('2022-01-03')
# current date = date that data were compiled

# 2 weeks to use as end dates for rolling 4-week bins
state_time_end = (current_day - as.numeric(format(current_day, '%w'))) - 15 - 0:1*(7)

# get the week number that corresponds to each date defined in state_time_end
data_weeks <- as.numeric(state_time_end - week0day1) %/% 7

# all combinations of states, 4-wk periods, & variants
# (reverse column order for convenience)
all.state = expand.grid(Variant = voc,
                        Roll_4wk_end = data_weeks,
                        State = sort(unique(dat$STUSAB)))[, 3:1]

# calculate estimated proportions (and CI) using survey design
ests = apply(X = all.state,
             MARGIN = 1,
             FUN = function(rr) myciprop(voc = rr[3],
                                         geoid = rr[1],
                                         svy = subset(svyDES,
                                                      week >= (as.numeric(rr[2]) - 3) &
                                                        week < (as.numeric(rr[2]) + 1)),
                                         str = FALSE))

# add together the combinations of states, time period, & variants with estimates
all.state = cbind(all.state,
                  Share    = ests[1,],
                  Share_lo = ests[2,],
                  Share_hi = ests[3,],
                  DF       = ests[4,],
                  eff.size = ests[5,],
                  cv.mean  = ests[6,])

# all combinations of states, 4-wk periods for "Other" variants
# (column order reversed for convenience)
others = expand.grid(Variant = "Other",
                     Roll_4wk_end = data_weeks,
                     State = sort(unique(dat$STUSAB)))[, 3:1]

# calculate the estimated proportions (and CI) using survey design
ests.others = apply(X = others,
                    MARGIN = 1,
                    FUN = function(rr) myciprop(voc = voc,
                                                geoid = rr[1],
                                                svy = subset(svyDES,
                                                             week >= (as.numeric(rr[2]) - 3) &
                                                               week < (as.numeric(rr[2]) + 1)),
                                                str = FALSE))

# add together the combinations of states, time period, & variants with estimates
others = cbind(others,
               Share    = 1-ests.others[1,],
               Share_lo = 1-ests.others[3,],
               Share_hi = 1-ests.others[2,],
               DF       = ests.others[4,],
               eff.size = ests.others[5,],
               cv.mean  = ests.others[6,])

# combine estimates for individual variants with estimates for "Other" variants
all.state = rbind(all.state,
                  others)

# empty object to hold sequence counts
all.state.out <-c()

# generate sequence counts by lineage, location, & date
for(i in seq(data_weeks)){
  
  # subset the data to the relevant time period
  dat2 <- dat[dat$week >= (as.numeric(data_weeks[i]) - 3) &
                dat$week < (as.numeric(data_weeks[i]) + 1),]
  
  # raw sample counts by variant & state
  raw_counts_state <- aggregate(count ~ VARIANT2 + STUSAB,
                                data = dat2,
                                FUN  = sum,
                                drop = FALSE)
  
  # set NA counts to 0
  raw_counts_state$count <- ifelse(test = is.na(raw_counts_state$count) == T,
                                   yes = 0,
                                   no = raw_counts_state$count)
  
  
  # merge sequence counts with weighted proportions estimates
  all.state2 <- merge(x = all.state[all.state$Roll_4wk_end == data_weeks[i],],
                      y = raw_counts_state,
                      by.x = c("State",
                               "Variant"),
                      by.y = c("STUSAB",
                               "VARIANT2"),
                      all = T)
  
  # set NA counts to 0 again
  all.state2[is.na(all.state2$count) == T, "count"] <- 0
  
  # calculate denominator counts
  dss <- aggregate(count ~ State,
                   data = all.state2,
                   FUN = sum)
  
  # change the column name for the denominator counts
  names(dss)[grep("count",names(dss))] <- "denom_count"
  
  # add the denominator counts into the results dataframe
  all.state2 <- merge(x = all.state2,
                      y = dss,
                      all = T)
  
  # add a column for the rolling 4-week period
  all.state2$Roll_Fourweek_ending <- unique(as.Date(dat$yr_wk[dat$week==data_weeks[i]]) + 6)
  
  # add the results for the given time period onto the dataframe of results
  all.state.out <- rbind.data.frame(all.state.out,
                                    all.state2)
} # end loop over weeks

# set the proportion (i.e. Share) to 0 and CI limits to NA when the count for a lineage is 0
all.state.out$Share = ifelse(test = all.state.out$Share != 0 &
                               all.state.out$count == 0,
                             yes = 0,
                             no = all.state.out$Share)
all.state.out$Share_lo = ifelse(test = is.na(all.state.out$Share_lo) == F &
                                  all.state.out$count == 0,
                                yes = NA,
                                no = all.state.out$Share_lo)
all.state.out$Share_hi = ifelse(test = is.na(all.state.out$Share_hi) == F &
                                  all.state.out$count == 0,
                                yes = NA,
                                no = all.state.out$Share_hi)

#calculate absolute CI width
all.state.out$CI_width = all.state.out$Share_hi - all.state.out$Share_lo

## generate NCHS flags
# flag estimates with Degrees of Freedom < 8
all.state.out$flag_df = ifelse(all.state.out$DF<8, 1, 0)
# flag estimates with effective size of < 30 (or NA)
all.state.out$flag_eff.size = ifelse(test = all.state.out$eff.size < 30 |
                                       is.na(all.state.out$eff.size)==T,
                                     yes = 1,
                                     no = 0)
# flag estimates with "denominator count" of < 30 (or NA)
# (denominator = count of sequences of all variants in a given region and time period)
all.state.out$flag_dss = ifelse(test = all.state.out$denom_count < 30|
                                  is.na(all.state.out$denom_count) == T,
                                yes = 1,
                                no = 0)
# flag estimates with wide (absolute) confidence intervals
all.state.out$flag_abs.ciw = ifelse(test = all.state.out$CI_width > 0.30 |
                                      is.na(all.state.out$CI_width) == T,
                                    yes = 1,
                                    no = 0)
# flag estimates with wide (relative) confidence intervals
all.state.out$flag_rel.ciw = ifelse(test = ((all.state.out$CI_width/all.state.out$Share)*100) > 130 |
                                      is.na((all.state.out$CI_width/all.state.out$Share)*100) == T,
                                    yes = 1,
                                    no = 0)

# Single identifier for observations that have *any* NCHS flag
all.state.out$nchs_flag = ifelse(test = all.state.out$flag_df == 1 |
                                   all.state.out$flag_eff.size == 1 |
                                   all.state.out$denom_count == 1 |
                                   all.state.out$flag_abs.ciw == 1 |
                                   all.state.out$flag_rel.ciw==1,
                                 yes = 1,
                                 no = 0)

# Single identifier for observations that have any NCHS flag *other* than the degrees of freedom flag.
all.state.out$nchs_flag_wodf = ifelse(test = all.state.out$flag_eff.size == 1 |
                                        all.state.out$denom_count == 1 |
                                        all.state.out$flag_abs.ciw == 1 |
                                        all.state.out$flag_rel.ciw == 1,
                                      yes = 1,
                                      no = 0)

# select columns for the final results
all.state.out = all.state.out[,c("State",
                                 "Roll_Fourweek_ending",
                                 "Variant",
                                 "Share",
                                 "Share_lo",
                                 "Share_hi",
                                 "count",
                                 "denom_count",
                                 "DF",
                                 "eff.size",
                                 "CI_width",
                                 "nchs_flag",
                                 "nchs_flag_wodf")]

# write results to file
write.csv(x = all.state.out,
          file = "state_weighted_roll4wk.csv",
          row.names = FALSE)
