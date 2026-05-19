# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Modified Variant Surveillance Code
#
# Adapted from MMWR Genomic Surveillance code (Paul, 2021-05-11 / 2022-01-18)
#
# Modifications:
#   1. Estimate period changed from 2 weeks (fortnight) to 4 weeks.
#   2. Weighted proportion CIs now use prop.test on weighted counts instead
#      of the survey-design-based svycipropkg (Korn-Graubard).
#   3. All estimates are calculated for the USA only; HHS region estimates
#      and the regional nowcast model have been removed.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Part 1: Define Functions ----------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## prop.test-based weighted proportion CI -------------------------------------
# Replacement for the survey-design-based myciprop/svycipropkg.
#
# This computes a weighted point estimate and a binomial confidence interval
# using stats::prop.test on the weighted counts. The effective "trials" is the
# sum of survey weights for the subset of interest, and the effective
# "successes" is the sum of weights for sequences belonging to the focal
# variant(s). The CI therefore reflects binomial uncertainty around the
# weighted proportion rather than full survey-design variance.
#
# Arguments:
#   ~ voc:   character vector of variant(s) to group together as the numerator
#   ~ dat:   data.table or data.frame containing VARIANT and weights columns,
#            already subset to the time period of interest (or pass full data
#            and use the `subset_expr` argument)
#   ~ weight_col: name of the column holding survey weights (default "weights")
#   ~ variant_col: name of the variant column (default "VARIANT")
#   ~ level: confidence level (default 0.95)
#   ~ mut: if TRUE, match against S_MUT mutation profile instead of VARIANT
#
# Output: named numeric vector with:
#   estimate, lcl, ucl, n_seq (raw sequence count in denominator),
#   n_voc (raw sequence count in numerator), weighted_denom, weighted_num
proptest_ci <- function(voc,
                        dat,
                        weight_col = "weights",
                        variant_col = "VARIANT",
                        level = 0.95,
                        mut = FALSE) {

  # coerce to data.frame to make column access uniform
  dat <- as.data.frame(dat)

  # if there are no rows, return NAs
  if (nrow(dat) == 0) {
    return(c(estimate = NA_real_,
             lcl = NA_real_,
             ucl = NA_real_,
             n_seq = 0,
             n_voc = 0,
             weighted_denom = 0,
             weighted_num = 0))
  }

  # binary indicator for membership in the focal variant group
  if (mut) {
    VOC <- grepl(
      pattern = paste0("^(?=.*\\",
                       paste(voc, collapse = "\\b)(?=.*\\"), "\\b)"),
      x = dat$S_MUT,
      perl = TRUE
    )
  } else {
    VOC <- (dat[[variant_col]] %in% voc)
  }

  # pull weights
  w <- dat[[weight_col]]

  # raw sequence counts
  n_seq <- nrow(dat)
  n_voc <- sum(VOC)

  # weighted counts
  weighted_denom <- sum(w, na.rm = TRUE)
  weighted_num   <- sum(w[VOC], na.rm = TRUE)

  # if there's no denominator weight or no sequences, return NAs for the CI
  if (weighted_denom <= 0 || n_seq == 0) {
    return(c(estimate = NA_real_,
             lcl = NA_real_,
             ucl = NA_real_,
             n_seq = n_seq,
             n_voc = n_voc,
             weighted_denom = weighted_denom,
             weighted_num = weighted_num))
  }

  # point estimate of the weighted proportion
  p_hat <- weighted_num / weighted_denom

  # prop.test requires integer x and n. Round weighted counts to nearest
  # integer for the CI calculation. The point estimate is preserved from
  # the un-rounded weighted proportion above; prop.test is used only for
  # the CI bounds.
  x <- round(weighted_num)
  n <- round(weighted_denom)

  # guard against x > n after rounding (rare edge case)
  if (x > n) x <- n

  # edge cases that prop.test does not handle gracefully
  if (n == 0) {
    lcl <- NA_real_
    ucl <- NA_real_
  } else if (p_hat == 0) {
    # use prop.test with x = 0 to get the one-sided upper bound
    pt <- suppressWarnings(prop.test(x = 0, n = n, conf.level = level))
    lcl <- 0
    ucl <- pt$conf.int[2]
  } else if (p_hat == 1) {
    pt <- suppressWarnings(prop.test(x = n, n = n, conf.level = level))
    lcl <- pt$conf.int[1]
    ucl <- 1
  } else {
    pt <- suppressWarnings(prop.test(x = x, n = n, conf.level = level))
    lcl <- pt$conf.int[1]
    ucl <- pt$conf.int[2]
  }

  c(estimate = p_hat,
    lcl = lcl,
    ucl = ucl,
    n_seq = n_seq,
    n_voc = n_voc,
    weighted_denom = weighted_denom,
    weighted_num = weighted_num)
}


## ---------------------------------------------------------------------------
## Nowcast functions (svymultinom, se.multinom, svyCI) -- unchanged from the
## original code. Bring these in verbatim from variant_surveillance_system.R:
##   - svymultinom
##   - se.multinom
##   - svyCI
##   - np / nearest_parent
##   - `%notin%`
##
## They're omitted here for brevity. Source the original file or copy them
## directly. Only the svycipropkg / myciprop functions are no longer needed
## for the weighted estimates, though svyCI is still used inside the nowcast
## pipeline.
## ---------------------------------------------------------------------------

source("variant_surveillance_system.R", local = FALSE,
       # If you don't want to re-run the original analysis, you can instead
       # copy just the function definitions (Part 1) into this file.
       echo = FALSE)
# NOTE: sourcing the original file will also execute Parts 2 and 3 of the
# original analysis. If that is not desired, copy just the function block
# (lines ~1-913 of variant_surveillance_system.R) into this file in place
# of the source() call above.


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Part 2: Data Prep & weight calculation --------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(survey)
library(nnet)
library(data.table)

options(survey.adjust.domain.lonely = TRUE,
        survey.lonely.psu = "average",
        stringsAsFactors = FALSE)

## Load sequence data ---------------------------------------------------------
dat <- readRDS(file = "../example/example_data.RDS")
dat$count <- 1
dat <- data.table::as.data.table(dat)

## Parameters -----------------------------------------------------------------
week0day1 <- as.Date("2020-01-05")

voc <- c("BA.1.1", "BA.2", "BA.2.12.1", "BA.2.75", "BA.2.75.2",
         "CH.1.1", "BF.7", "BF.11", "BA.4", "BA.4.6", "BA.5", "BA.5.2.6",
         "BQ.1", "BQ.1.1", "BN.1", "XBB", "XBB.1.5", "XBB.1.5.1",
         "XBB.1.9.1", "FD.2", "XBB.1.9.2", "XBB.1.16", "XBB.2.3",
         "B.1.617.2", "B.1.1.529")

## Aggregate sublineages ------------------------------------------------------
dat$VARIANT <- dat$lineage <- as.character(dat$VARIANT)
lut <- dat[, .(expanded_lineage = unique(expanded_lineage)), by = "VARIANT"]
voc_expanded <- lut[VARIANT %in% voc][["expanded_lineage"]]
unique_vars <- na.omit(unique(dat$VARIANT))
unique_lineage_expanded <- unname(setNames(lut$expanded_lineage,
                                           lut$VARIANT)[unique_vars])

voc_lut <- data.frame(
  variant = unique_vars,
  lineage_expanded = unique_lineage_expanded,
  parent_lineage_expanded = nearest_parent(unique_lineage_expanded,
                                           voc_expanded)
)
voc_lut$parent_variant <- setNames(voc_lut$variant,
                                   voc_lut$lineage_expanded)[voc_lut$parent_lineage_expanded]
row.names(voc_lut) <- 1:nrow(voc_lut)

dat$VARIANT <- setNames(voc_lut$parent_variant, voc_lut$variant)[dat$VARIANT]
dat$VARIANT2 <- as.character(dat$VARIANT)
dat[dat$VARIANT %notin% voc, "VARIANT2"] <- "Other"

## Survey weights (unchanged) -------------------------------------------------
dat[, "proxy_infections" :=
      state_population *
      sqrt(POSITIVE.HHS.nrevss / TOTAL.HHS.nrevss) *
      sqrt(POSITIVE.HHS / population_reporting.HHS)]

dat[, "weight" := proxy_infections / sum(count),
    by = c("STUSAB", "yr_wk")]

invalid_weight <- is.na(dat$weight) | is.infinite(dat$weight)
dat <- subset(dat, !invalid_weight)

## --- CHANGE 1: 4-week periods ------------------------------------------------
# The input data is assumed to already contain a FOURWEEK_END column that
# tags each sequence with the end date of its four-week estimate period.
# We simply coerce it to Date class for safe downstream date arithmetic.
if (!"FOURWEEK_END" %in% names(dat)) {
  stop("Input data must contain a FOURWEEK_END column.")
}
dat$FOURWEEK_END <- as.Date(dat$FOURWEEK_END)

## Survey design & weight trimming --------------------------------------------
data_date <- as.Date("2023-05-13")
current_week <- as.numeric(data_date - week0day1) %/% 7
dat$current_week <- current_week

svyDES <- survey::svydesign(ids = ~ SOURCE,
                            strata  = ~ STUSAB + yr_wk,
                            weights = ~ weight,
                            nest = TRUE,
                            data = dat)

max_weight <- quantile(weights(svyDES), probs = 0.99)
wts <- weights(svyDES)
min_wt <- min(wts[wts > 0])
svyDES <- survey::trimWeights(svyDES, upper = max_weight, lower = min_wt)
dat$weights <- weights(svyDES)

## Model-week centering (unchanged) -------------------------------------------
time_end <- "2023-05-13"
model_weeks <- 21
model_week_max <- as.numeric(as.Date(time_end) - week0day1) %/% 7
model_week_min <- model_week_max - model_weeks
model_week_mid <- round(model_weeks / 2)
dat$model_week <- dat$week - model_week_min - model_week_mid


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Part 3: Calculating Variant Proportions -------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## --- CHANGE 2 + 3: USA-only, prop.test-based weighted proportions ----------

# Get the four-week periods to estimate over
fwks <- sort(unique(dat$FOURWEEK_END))

# Combinations: only USA, one row per variant per 4-week period
all.fwk <- expand.grid(
  Variant = voc,
  Fourweek_ending = fwks,
  USA_or_HHSRegion = "USA",
  stringsAsFactors = FALSE
)[, 3:1]

# Convenience: keep weights on the data.table for prop.test computations
dat_for_props <- as.data.frame(dat)

# Estimate weighted proportion + prop.test CI for each variant in each period
ests <- apply(
  X = all.fwk,
  MARGIN = 1,
  FUN = function(rr) {
    # rr[1] = USA_or_HHSRegion (always "USA")
    # rr[2] = Fourweek_ending (as character; coerce back to Date)
    # rr[3] = Variant
    sub <- dat_for_props[dat_for_props$FOURWEEK_END == as.Date(rr[2]), ]
    proptest_ci(voc = rr[3],
                dat = sub,
                level = 0.95)
  }
)

all.fwk <- cbind(
  all.fwk,
  Share          = ests["estimate", ],
  Share_lo       = ests["lcl", ],
  Share_hi       = ests["ucl", ],
  n_seq          = ests["n_seq", ],
  n_voc          = ests["n_voc", ],
  weighted_denom = ests["weighted_denom", ],
  weighted_num   = ests["weighted_num", ]
)

# "Other" estimates: compute the aggregate proportion across all voc, then
# subtract from 1 to get the residual "Other" proportion and flip the CI.
others <- expand.grid(
  Variant = "Other",
  Fourweek_ending = fwks,
  USA_or_HHSRegion = "USA",
  stringsAsFactors = FALSE
)[, 3:1]

ests_others <- apply(
  X = others,
  MARGIN = 1,
  FUN = function(rr) {
    sub <- dat_for_props[dat_for_props$FOURWEEK_END == as.Date(rr[2]), ]
    proptest_ci(voc = voc, dat = sub, level = 0.95)
  }
)

others <- cbind(
  others,
  Share          = 1 - ests_others["estimate", ],
  Share_lo       = 1 - ests_others["ucl", ],
  Share_hi       = 1 - ests_others["lcl", ],
  n_seq          = ests_others["n_seq", ],
  n_voc          = ests_others["n_seq", ] - ests_others["n_voc", ],
  weighted_denom = ests_others["weighted_denom", ],
  weighted_num   = ests_others["weighted_denom", ] - ests_others["weighted_num", ]
)

all.fwk2 <- rbind(all.fwk, others)

# Add CI width and NCHS-style reliability flags. With prop.test there is no
# survey degrees-of-freedom or design effect to flag on, but we can still
# flag small denominators and wide CIs.
all.fwk2$CI_width      <- all.fwk2$Share_hi - all.fwk2$Share_lo
all.fwk2$flag_n_seq    <- as.numeric(all.fwk2$n_seq < 30 | is.na(all.fwk2$n_seq))
all.fwk2$flag_abs_ciw  <- as.numeric(all.fwk2$CI_width > 0.30 |
                                       is.na(all.fwk2$CI_width))
all.fwk2$flag_rel_ciw  <- as.numeric(
  (all.fwk2$CI_width / all.fwk2$Share) * 100 > 130 |
    is.na((all.fwk2$CI_width / all.fwk2$Share) * 100)
)
all.fwk2$any_flag <- as.numeric(
  all.fwk2$flag_n_seq == 1 |
    all.fwk2$flag_abs_ciw == 1 |
    all.fwk2$flag_rel_ciw == 1
)

write.csv(all.fwk2,
          file = "fourweekly_weighted_proportions.csv",
          row.names = FALSE)


## --- Nowcast: USA-only --------------------------------------------------------

# Variants to include in the nowcast model (same logic as the original)
n_recent_weeks <- 4
us_var <- sort(
  prop.table(xtabs(weights ~ VARIANT,
                   data = dat,
                   subset = (dat$week >= current_week - n_recent_weeks))),
  decreasing = TRUE
)
model_vars <- voc[voc %in% dat[model_week %in% ((1:model_weeks) - model_week_mid)][["VARIANT"]]]

dat$K_US <- match(dat$VARIANT, model_vars)
dat$K_US[is.na(dat$K_US)] <- length(model_vars) + 1

moddat <- subset(dat, model_week %in% ((1:model_weeks) - model_week_mid))
moddat$wts <- moddat$weights / max(moddat$weights)

mysvy <- survey::svydesign(ids = ~ SOURCE,
                           strata  = ~ STUSAB + yr_wk,
                           weights = ~ wts,
                           nest = TRUE,
                           data = moddat)

# Fit USA-only nowcast model (no HHS predictor)
svymlm_us <- svymultinom(
  mod.dat    = moddat,
  mysvy      = mysvy,
  fmla       = formula("as.numeric(as.factor(K_US)) ~ model_week"),
  model_vars = model_vars
)

## Aggregation matrix (unchanged) ---------------------------------------------
output_lineages <- model_vars[model_vars %in% voc]
agg_var_mat <- matrix(0, nrow = 0, ncol = (length(model_vars) + 1))
colnames(agg_var_mat) <- c(model_vars, "Other")

model_vars_expanded <- setNames(voc_lut$lineage_expanded,
                                voc_lut$variant)[model_vars]
output_lineages_expanded <- setNames(voc_lut$lineage_expanded,
                                     voc_lut$variant)[output_lineages]
model_var_parents_expanded <- nearest_parent(model_vars_expanded,
                                             output_lineages_expanded)
model_var_parents <- setNames(voc_lut$variant,
                              voc_lut$lineage_expanded)[model_var_parents_expanded]

model_var_lut <- data.frame(
  model_vars = model_vars,
  output_model_vars = model_var_parents
)
mvl_sub <- model_var_lut[model_var_lut$model_vars != model_var_lut$output_model_vars, ]
unique_mvl_sub <- unique(mvl_sub$output_model_vars)

for (i in unique_mvl_sub) {
  extra_row <- ifelse(colnames(agg_var_mat) %in%
                        c(i, model_var_lut[model_var_lut$output_model_vars == i, "model_vars"]),
                      1, 0)
  agg_var_mat <- rbind(agg_var_mat, extra_row)
  row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(i, "Aggregated")
}
other_agg <- base::setdiff(colnames(agg_var_mat)[colSums(agg_var_mat) == 0],
                           output_lineages)
agg_var_mat <- rbind(
  agg_var_mat,
  ifelse(colnames(agg_var_mat) %in% other_agg, 1, 0)
)
row.names(agg_var_mat)[nrow(agg_var_mat)] <- "Other Aggregated"


## Four-weekly nowcast projections (USA only) ---------------------------------
# Use the final 4-week period of observed data plus two 4-week periods ahead.
proj_fwks <- as.Date(tail(fwks, 1))
proj_fwks <- sort(unique(c(proj_fwks,
                           proj_fwks + 28,
                           proj_fwks + 56)))

proj.res <- c()

for (fwk in proj_fwks) {

  mlm <- svymlm_us$mlm
  geoid <- "USA"

  # Use the midpoint of the four-week period as the prediction time.
  # A 4-week window has its midpoint 13.5 days before its end date.
  wk <- as.numeric(as.Date(fwk, origin = "1970-01-01") -
                     (week0day1 + 3) - 13.5) / 7
  wk <- wk - model_week_min - model_week_mid

  ests <- se.multinom(mlm = mlm,
                      newdata_1row = data.frame(model_week = wk),
                      composite_variant = agg_var_mat)

  gr <- with(ests, 100 * exp(b_i - sum(p_i * b_i)) - 100)
  se.gr <- with(ests,
                100 * exp(sqrt(se.b_i^2 * (1 - 2 * p_i) +
                                 sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)

  se.gr_link <- with(ests,
                     sqrt(se.b_i^2 * (1 - 2 * p_i) +
                            sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))
  gr_link <- with(ests, (b_i - sum(p_i * b_i)))
  gr_lo_link <- gr_link - 1.96 * se.gr_link
  gr_hi_link <- gr_link + 1.96 * se.gr_link

  gr_lo <- 100 * exp(gr_lo_link) - 100
  gr_hi <- 100 * exp(gr_hi_link) - 100

  doubling_time    <- log(2) / gr_link    * 7
  doubling_time_lo <- log(2) / gr_lo_link * 7
  doubling_time_hi <- log(2) / gr_hi_link * 7

  gr_agg <- data.frame(variant = rownames(agg_var_mat),
                       gr = NA, se.gr = NA, gr_lo = NA, gr_hi = NA,
                       dt = NA, dt_lo = NA, dt_hi = NA)
  for (r in 1:nrow(agg_var_mat)) {
    if (unname(rowSums(agg_var_mat)[r]) == 1) {
      col_ind <- which(agg_var_mat[r, ] > 0)
      gr_agg[r, "gr"]    <- gr[col_ind]
      gr_agg[r, "se.gr"] <- se.gr[col_ind]
      gr_agg[r, "gr_lo"] <- gr_lo[col_ind]
      gr_agg[r, "gr_hi"] <- gr_hi[col_ind]
      gr_agg[r, "dt"]    <- doubling_time[col_ind]
      gr_agg[r, "dt_lo"] <- doubling_time_lo[col_ind]
      gr_agg[r, "dt_hi"] <- doubling_time_hi[col_ind]
    } else {
      col_ind <- unname(which(agg_var_mat[r, ] > 0))
      gr_agg[r, "gr"]    <- sum(gr[col_ind]    * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
      gr_agg[r, "gr_lo"] <- sum(gr_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
      gr_agg[r, "gr_hi"] <- sum(gr_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
      gr_agg[r, "dt"]    <- sum(doubling_time[col_ind]    * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
      gr_agg[r, "dt_lo"] <- sum(doubling_time_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
      gr_agg[r, "dt_hi"] <- sum(doubling_time_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
    }
  }

  ests_dt <- data.table::data.table(
    USA_or_HHSRegion = geoid,
    Fourweek_ending  = as.Date(fwk, origin = "1970-01-01"),
    Variant = c(model_vars, "Other",
                row.names(ests$composite_variant$matrix)),
    Share    = c(ests$p_i, ests$composite_variant$p_i),
    se.Share = c(ests$se.p_i, ests$composite_variant$se.p_i),
    growth_rate    = c(gr,    gr_agg$gr),
    growth_rate_lo = c(gr_lo, gr_agg$gr_lo),
    growth_rate_hi = c(gr_hi, gr_agg$gr_hi),
    doubling_time    = c(doubling_time,    gr_agg$dt),
    doubling_time_lo = c(doubling_time_lo, gr_agg$dt_lo),
    doubling_time_hi = c(doubling_time_hi, gr_agg$dt_hi)
  )

  binom.ci <- apply(ests_dt, 1,
                    function(rr) svyCI(p = as.numeric(rr[4]),
                                       s = as.numeric(rr[5]),
                                       conf.level = 0.95))
  ests_dt$Share_lo <- binom.ci[1, ]
  ests_dt$Share_hi <- binom.ci[2, ]

  proj.res <- rbind(proj.res, ests_dt)
}

## Write nowcast outputs ------------------------------------------------------
agg_lineages <- colnames(agg_var_mat)[colSums(agg_var_mat) > 0]
if ("Other" %notin% agg_lineages) agg_lineages <- c(agg_lineages, "Other Aggregated")
results_agg <- proj.res[Variant %notin% agg_lineages]

if (!all(results_agg[, .(total_share = sum(Share)),
                     by = c("USA_or_HHSRegion", "Fourweek_ending")
                     ][, unique(round(total_share, 5))] == 1)) {
  warning("Total proportion does not add up to 100% in each time period (aggregated).")
} else {
  write.csv(results_agg,
            file = "fourweekly_nowcast_proportions_aggregated.csv",
            row.names = FALSE)
}

drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]
results_nonagg <- proj.res[
  Variant %notin% c(drop_lin, "Other",
                    colnames(agg_var_mat["Other Aggregated", , drop = FALSE])[
                      agg_var_mat["Other Aggregated", , drop = FALSE] > 0])
]

if (!all(results_nonagg[, .(total_share = sum(Share)),
                        by = c("USA_or_HHSRegion", "Fourweek_ending")
                        ][, unique(round(total_share, 5))] == 1)) {
  warning("Total proportion does not add up to 100% in each time period (non-aggregated).")
} else {
  write.csv(results_nonagg,
            file = "fourweekly_nowcast_proportions.csv",
            row.names = FALSE)
}