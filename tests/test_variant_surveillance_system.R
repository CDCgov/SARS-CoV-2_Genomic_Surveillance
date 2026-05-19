# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Unit tests for variant_surveillance_modified.R
#
# Run with:
#   testthat::test_file("test_variant_surveillance.R")
#
# Assumptions:
#   - variant_surveillance_system.R is on the working directory (sourced for
#     the unchanged helper functions: svymultinom, se.multinom, svyCI,
#     nearest_parent, np, %notin%).
#   - variant_surveillance_modified.R defines proptest_ci.
#   - example_data.RDS is on the working directory and is used for the
#     integration-style tests of svymultinom / se.multinom.
#
# These tests are organized in three layers:
#   1. Pure-function unit tests on synthetic inputs (fast, deterministic).
#   2. proptest_ci tests with a small hand-built data.frame.
#   3. Integration tests on example_data.RDS that check invariants of the
#      multinomial nowcast pipeline (proportions sum to 1, SEs are
#      non-negative, dimensions are right). These are slower.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(testthat)
library(survey)
library(nnet)
library(data.table)

# Source the helper functions. Wrap in tryCatch so the test file gives a clear
# error if the source files are missing rather than failing inside the first
# test_that() block.
if (!file.exists("variant_surveillance_system.R")) {
  stop("variant_surveillance_system.R not found in working directory.")
}
if (!file.exists("variant_surveillance_modified.R")) {
  stop("variant_surveillance_modified.R not found in working directory.")
}

# To avoid running Parts 2 and 3 of the original script during sourcing,
# extract just the function definitions. The simplest approach: source the
# original file in a temporary environment but suppress side effects by
# overriding readRDS to return a tiny placeholder. If that's not workable
# in your setup, copy the function block into a separate helpers.R file.
#
# Here we use a more direct approach: source the file but immediately catch
# any error from the data-loading step. Functions defined before the error
# will still be available.
.helper_env <- new.env()
tryCatch(
  sys.source("variant_surveillance_system.R", envir = .helper_env),
  error = function(e) {
    message("Note: original script halted partway through sourcing (expected ",
            "if example_data.RDS or downstream objects are unavailable). ",
            "Continuing with whatever functions were successfully defined.")
  }
)

# Pull the helpers we need into the global env for the tests
for (fn in c("svymultinom", "se.multinom", "svyCI", "nearest_parent", "np",
             "%notin%", "myciprop", "svycipropkg")) {
  if (exists(fn, envir = .helper_env, inherits = FALSE)) {
    assign(fn, get(fn, envir = .helper_env), envir = .GlobalEnv)
  }
}

# Source the modified script's function block. proptest_ci is the only new
# function; everything else is unchanged. We avoid running its Parts 2-3 by
# extracting just the proptest_ci definition.
.modified_lines <- readLines("variant_surveillance_modified.R")
.start <- grep("^proptest_ci\\s*<-\\s*function", .modified_lines)
if (length(.start) == 0) {
  stop("Could not locate proptest_ci definition in variant_surveillance_modified.R")
}
# Find the matching closing brace at column 1
.brace_close <- grep("^\\}\\s*$", .modified_lines)
.end <- min(.brace_close[.brace_close > .start])
eval(parse(text = .modified_lines[.start:.end]), envir = .GlobalEnv)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Layer 1: pure-function unit tests on synthetic inputs ----------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

test_that("%notin% returns the logical negation of %in%", {
  expect_equal(1 %notin% c(2, 3, 4), TRUE)
  expect_equal(2 %notin% c(2, 3, 4), FALSE)
  expect_equal(c("a", "b", "z") %notin% c("a", "b", "c"),
               c(FALSE, FALSE, TRUE))
  # empty haystack
  expect_equal("x" %notin% character(0), TRUE)
  # NA handling: NA %in% anything is NA, so %notin% should also be NA
  expect_equal(NA %notin% c(1, 2), NA)
})

test_that("np finds the longest matching parent lineage", {
  # BA.5.2.6 should resolve to BA.5 (longest parent in the candidate set)
  expect_equal(unname(np("BA.5.2.6", c("BA.5", "BA.2", "BA"))), "BA.5")
  # BQ.1.1 should resolve to BQ.1
  expect_equal(unname(np("BQ.1.1", c("BQ.1", "BQ"))), "BQ.1")
  # an exact match returns itself
  expect_equal(unname(np("XBB.1.5", c("XBB.1.5", "XBB"))), "XBB.1.5")
  # no match returns the no_match sentinel
  expect_equal(unname(np("AY.4", c("BA.1", "BA.2"))), "Other")
  expect_equal(unname(np("AY.4", c("BA.1"), no_match = "Unassigned")),
               "Unassigned")
})

test_that("np does not confuse partial prefix matches", {
  # BA.5.22 must not match BA.5.2 just because of leading characters;
  # the trailing-dot mechanism inside np() should prevent this.
  expect_equal(unname(np("BA.5.22", c("BA.5.2", "BA.5"))), "BA.5")
})

test_that("nearest_parent vectorizes np across multiple inputs", {
  res <- nearest_parent(c("BA.5.2.6", "BQ.1.1", "AY.4"),
                        c("BA.5", "BQ.1", "BA.2"))
  expect_equal(unname(res), c("BA.5", "BQ.1", "Other"))
  expect_equal(names(res), c("BA.5.2.6", "BQ.1.1", "AY.4"))
})

test_that("svyCI returns degenerate intervals at the boundaries", {
  # SE of zero -> CI of (0, 0) by convention in the function
  expect_equal(svyCI(p = 0.5, s = 0), c(0, 0))
  # p == 0 -> CI of (0, 0)
  expect_equal(svyCI(p = 0,   s = 0.01), c(0, 0))
  # p == 1 -> CI of (1, 1)
  expect_equal(svyCI(p = 1,   s = 0.01), c(1, 1))
})

test_that("svyCI returns a valid binomial CI in the interior", {
  ci <- svyCI(p = 0.3, s = 0.05, conf.level = 0.95)
  expect_length(ci, 2)
  expect_true(ci[1] >= 0 && ci[2] <= 1)
  expect_true(ci[1] < 0.3 && ci[2] > 0.3)
  # narrower CIs for smaller SE
  ci_narrow <- svyCI(p = 0.3, s = 0.01, conf.level = 0.95)
  expect_lt(ci_narrow[2] - ci_narrow[1], ci[2] - ci[1])
})


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Layer 2: proptest_ci on a small hand-built dataset -------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

make_toy <- function() {
  # 10 sequences, two variants, equal weights -> simple proportion of 0.4
  data.frame(
    VARIANT = c("BA.5", "BA.5", "BA.5", "BA.5",
                "XBB",  "XBB",  "XBB",  "XBB",  "XBB",  "XBB"),
    weights = rep(1, 10),
    stringsAsFactors = FALSE
  )
}

test_that("proptest_ci recovers the unweighted proportion when weights are equal", {
  d <- make_toy()
  out <- proptest_ci(voc = "BA.5", dat = d)
  expect_equal(unname(out["estimate"]), 0.4)
  expect_equal(unname(out["n_seq"]), 10)
  expect_equal(unname(out["n_voc"]), 4)
  expect_equal(unname(out["weighted_denom"]), 10)
  expect_equal(unname(out["weighted_num"]), 4)
})

test_that("proptest_ci handles unequal weights correctly", {
  d <- make_toy()
  # double-weight the BA.5 sequences
  d$weights <- ifelse(d$VARIANT == "BA.5", 2, 1)
  out <- proptest_ci(voc = "BA.5", dat = d)
  # weighted: numerator = 4*2 = 8, denominator = 8 + 6 = 14
  expect_equal(unname(out["weighted_num"]),   8)
  expect_equal(unname(out["weighted_denom"]), 14)
  expect_equal(unname(out["estimate"]), 8 / 14)
})

test_that("proptest_ci returns a valid binomial CI bracketing the estimate", {
  d <- make_toy()
  out <- proptest_ci(voc = "BA.5", dat = d, level = 0.95)
  expect_true(out["lcl"] <= out["estimate"])
  expect_true(out["ucl"] >= out["estimate"])
  expect_true(out["lcl"] >= 0 && out["ucl"] <= 1)
})

test_that("proptest_ci groups multiple variants into a single numerator", {
  d <- make_toy()
  out <- proptest_ci(voc = c("BA.5", "XBB"), dat = d)
  expect_equal(unname(out["estimate"]), 1.0)
  expect_equal(unname(out["n_voc"]), 10)
})

test_that("proptest_ci handles the p = 0 boundary", {
  d <- make_toy()
  out <- proptest_ci(voc = "BQ.1", dat = d)  # not in data
  expect_equal(unname(out["estimate"]), 0)
  expect_equal(unname(out["lcl"]), 0)
  expect_true(out["ucl"] > 0)            # one-sided upper bound exists
  expect_true(out["ucl"] < 1)
})

test_that("proptest_ci handles the p = 1 boundary", {
  d <- make_toy()
  d$VARIANT <- "BA.5"
  out <- proptest_ci(voc = "BA.5", dat = d)
  expect_equal(unname(out["estimate"]), 1)
  expect_equal(unname(out["ucl"]), 1)
  expect_true(out["lcl"] < 1)
  expect_true(out["lcl"] > 0)
})

test_that("proptest_ci returns NA-filled output for empty data", {
  d <- make_toy()[0, ]
  out <- proptest_ci(voc = "BA.5", dat = d)
  expect_true(is.na(out["estimate"]))
  expect_true(is.na(out["lcl"]))
  expect_true(is.na(out["ucl"]))
  expect_equal(unname(out["n_seq"]), 0)
})

test_that("proptest_ci confidence level controls CI width", {
  d <- make_toy()
  out95 <- proptest_ci(voc = "BA.5", dat = d, level = 0.95)
  out99 <- proptest_ci(voc = "BA.5", dat = d, level = 0.99)
  # 99% CI should be at least as wide as 95%
  w95 <- out95["ucl"] - out95["lcl"]
  w99 <- out99["ucl"] - out99["lcl"]
  expect_gte(w99, w95)
})

test_that("proptest_ci with mut = TRUE matches against S_MUT", {
  d <- data.frame(
    VARIANT = rep("X", 4),
    S_MUT   = c("L452R+E484K", "L452R", "E484K", "wildtype"),
    weights = rep(1, 4),
    stringsAsFactors = FALSE
  )
  # voc = c("L452R", "E484K") with mut = TRUE matches strings containing BOTH
  out <- proptest_ci(voc = c("L452R", "E484K"), dat = d, mut = TRUE)
  expect_equal(unname(out["n_voc"]), 1)  # only the first row matches
  expect_equal(unname(out["estimate"]), 0.25)
})


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Layer 3: integration tests on example_data.RDS -----------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# These are slower and check invariants of svymultinom + se.multinom on real
# data rather than exact numeric outputs. They are skipped if the data file
# is unavailable.

skip_if_no_data <- function() {
  if (!file.exists("example_data.RDS")) {
    testthat::skip("example_data.RDS not available; skipping integration tests.")
  }
}

# Cache the prepared survey design + moddat so we don't redo it in each test
prepare_modeling_objects <- function() {
  dat <- readRDS("example_data.RDS")
  dat$count <- 1
  dat <- data.table::as.data.table(dat)
  dat$VARIANT <- as.character(dat$VARIANT)

  voc <- c("BA.5", "BQ.1", "BQ.1.1", "XBB", "XBB.1.5",
           "B.1.617.2", "B.1.1.529")

  # minimal prep: skip the sublineage aggregation for speed; treat VARIANT
  # as-is and bucket non-voc into "Other"
  dat[VARIANT %notin% voc, VARIANT := "Other"]

  # infection weights (same formula as the script)
  dat[, proxy_infections :=
        state_population *
        sqrt(POSITIVE.HHS.nrevss / TOTAL.HHS.nrevss) *
        sqrt(POSITIVE.HHS / population_reporting.HHS)]
  dat[, weight := proxy_infections / sum(count), by = c("STUSAB", "yr_wk")]
  dat <- subset(dat, !is.na(weight) & !is.infinite(weight))

  week0day1 <- as.Date("2020-01-05")
  data_date <- as.Date("2023-05-13")
  current_week <- as.numeric(data_date - week0day1) %/% 7
  model_weeks <- 21
  model_week_max <- as.numeric(data_date - week0day1) %/% 7
  model_week_min <- model_week_max - model_weeks
  model_week_mid <- round(model_weeks / 2)
  dat$model_week <- dat$week - model_week_min - model_week_mid

  svyDES <- survey::svydesign(ids = ~ SOURCE,
                              strata = ~ STUSAB + yr_wk,
                              weights = ~ weight,
                              nest = TRUE,
                              data = dat)
  wts <- weights(svyDES)
  svyDES <- survey::trimWeights(svyDES,
                                upper = quantile(wts, 0.99),
                                lower = min(wts[wts > 0]))
  dat$weights <- weights(svyDES)

  model_vars <- voc[voc %in% dat[model_week %in% ((1:model_weeks) - model_week_mid)][["VARIANT"]]]
  dat$K_US <- match(dat$VARIANT, model_vars)
  dat$K_US[is.na(dat$K_US)] <- length(model_vars) + 1

  moddat <- subset(dat, model_week %in% ((1:model_weeks) - model_week_mid))
  moddat$wts <- moddat$weights / max(moddat$weights)

  mysvy <- survey::svydesign(ids = ~ SOURCE,
                             strata = ~ STUSAB + yr_wk,
                             weights = ~ wts,
                             nest = TRUE,
                             data = moddat)

  list(dat = dat, moddat = moddat, mysvy = mysvy,
       model_vars = model_vars,
       model_week_min = model_week_min,
       model_week_mid = model_week_mid)
}

# memoize once per test run
.cached_prep <- NULL
get_prep <- function() {
  if (is.null(.cached_prep)) {
    .cached_prep <<- prepare_modeling_objects()
  }
  .cached_prep
}

test_that("svymultinom fits a USA-only model and returns expected structure", {
  skip_if_no_data()
  p <- get_prep()

  fit <- svymultinom(
    mod.dat = p$moddat,
    mysvy   = p$mysvy,
    fmla    = formula("as.numeric(as.factor(K_US)) ~ model_week"),
    model_vars = p$model_vars
  )

  expect_true(is.list(fit))
  expect_true(all(c("mlm", "estimates", "variants") %in% names(fit)))
  expect_s3_class(fit$mlm, "multinom")

  # if Hessian was invertible, we expect SE + sandwich + scores
  if (!is.null(fit$SE)) {
    expect_true(is.numeric(fit$SE))
    expect_true(all(fit$SE >= 0))
    expect_true(is.matrix(fit$sandwich))
    # sandwich should be square
    expect_equal(nrow(fit$sandwich), ncol(fit$sandwich))
    # and symmetric (up to numerical tolerance)
    expect_equal(fit$sandwich, t(fit$sandwich), tolerance = 1e-8)
    # diagonal entries (variances) should be non-negative
    expect_true(all(diag(fit$sandwich) >= -1e-10))
  }
})

test_that("se.multinom returns proportions that sum to 1", {
  skip_if_no_data()
  p <- get_prep()

  fit <- svymultinom(
    mod.dat = p$moddat,
    mysvy   = p$mysvy,
    fmla    = formula("as.numeric(as.factor(K_US)) ~ model_week"),
    model_vars = p$model_vars
  )

  # predict at a single model_week within range
  ests <- se.multinom(mlm = fit$mlm,
                      newdata_1row = data.frame(model_week = 0))

  expect_equal(sum(ests$p_i), 1, tolerance = 1e-10)
  expect_true(all(ests$p_i >= 0))
  expect_true(all(ests$p_i <= 1))
  # number of classes = model_vars + 1 ("Other")
  expect_equal(length(ests$p_i), length(p$model_vars) + 1)
})

test_that("se.multinom standard errors are non-negative and finite", {
  skip_if_no_data()
  p <- get_prep()

  fit <- svymultinom(
    mod.dat = p$moddat,
    mysvy   = p$mysvy,
    fmla    = formula("as.numeric(as.factor(K_US)) ~ model_week"),
    model_vars = p$model_vars
  )

  ests <- se.multinom(mlm = fit$mlm,
                      newdata_1row = data.frame(model_week = 0))

  expect_true(all(ests$se.p_i >= 0))
  expect_true(all(is.finite(ests$se.p_i)))
  expect_equal(length(ests$se.p_i), length(ests$p_i))
})

test_that("se.multinom composite_variant aggregates correctly", {
  skip_if_no_data()
  p <- get_prep()

  fit <- svymultinom(
    mod.dat = p$moddat,
    mysvy   = p$mysvy,
    fmla    = formula("as.numeric(as.factor(K_US)) ~ model_week"),
    model_vars = p$model_vars
  )

  # build a composite matrix: aggregate the first two model variants
  K <- length(p$model_vars) + 1
  cv <- matrix(0, nrow = 1, ncol = K)
  cv[1, 1:2] <- 1

  ests <- se.multinom(mlm = fit$mlm,
                      newdata_1row = data.frame(model_week = 0),
                      composite_variant = cv)

  # aggregated proportion equals sum of components
  expect_equal(unname(ests$composite_variant$p_i),
               sum(ests$p_i[1:2]),
               tolerance = 1e-10)
  # aggregated SE is non-negative
  expect_true(ests$composite_variant$se.p_i >= 0)
  # aggregated proportion in [0, 1]
  expect_true(ests$composite_variant$p_i >= 0 &&
              ests$composite_variant$p_i <= 1)
})

test_that("se.multinom predictions are stable across nearby weeks", {
  skip_if_no_data()
  p <- get_prep()

  fit <- svymultinom(
    mod.dat = p$moddat,
    mysvy   = p$mysvy,
    fmla    = formula("as.numeric(as.factor(K_US)) ~ model_week"),
    model_vars = p$model_vars
  )

  e0 <- se.multinom(fit$mlm, data.frame(model_week = 0))
  e_small <- se.multinom(fit$mlm, data.frame(model_week = 0.01))

  # small change in week -> small change in proportions
  expect_true(max(abs(e0$p_i - e_small$p_i)) < 0.05)
})
