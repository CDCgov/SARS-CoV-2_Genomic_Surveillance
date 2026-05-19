# Unit Tests

This folder contains the automated tests for the variant surveillance code. The
tests verify that the functions in `variant_surveillance_modified.R` (and the
unchanged helper functions inherited from `variant_surveillance_system.R`)
behave correctly.

## Requirements

- **R** (version 4.0 or later recommended)
- The following R packages:
  - `testthat` — the testing framework
  - `survey` — survey design functions used by the nowcast model
  - `nnet` — multinomial regression used by the nowcast model
  - `data.table` — used throughout the analysis code

Install the packages with:

```r
install.packages(c("testthat", "survey", "nnet", "data.table"))
```

### Files that must be available

The test script expects the following files to be reachable from the working
directory when the tests are run:

| File | Purpose | Required? |
|------|---------|-----------|
| `variant_surveillance_system.R` | Source of the unchanged helper functions (`svymultinom`, `se.multinom`, `svyCI`, `nearest_parent`, `np`, `%notin%`) | Yes |
| `variant_surveillance_modified.R` | Source of the new `proptest_ci` function | Yes |
| `example_data.RDS` | Real data used by the integration tests | Optional — see below |

If `example_data.RDS` is not present, the integration tests (Layer 3, described
below) are automatically skipped and the rest of the tests still run.

## What the tests cover

The tests are organized into three layers, from fastest and most isolated to
slowest and most integrated.

### Layer 1 — Pure-function unit tests

Fast, deterministic tests that run on small synthetic inputs and need no data
file. They cover the helper functions:

- **`%notin%`** — confirms it is the logical negation of `%in%`, including
  empty inputs and `NA` handling.
- **`np` / `nearest_parent`** — confirms each lineage resolves to its longest
  matching parent, that exact matches resolve to themselves, that the
  no-match sentinel is returned when appropriate, and that partial prefix
  matches are not mistaken for true parents.
- **`svyCI`** — confirms the boundary cases (zero standard error, proportion of
  0, proportion of 1) return the expected degenerate intervals, and that
  interior intervals are valid and narrow as the standard error shrinks.

### Layer 2 — `proptest_ci` tests

Tests for the new `proptest_ci` function, run against a small hand-built data
frame where the correct answers are known by construction:

- Recovers the unweighted proportion when all weights are equal.
- Computes the weighted numerator and denominator correctly when weights
  differ.
- Returns a binomial confidence interval that brackets the point estimate and
  stays within [0, 1].
- Groups multiple variants into a single numerator when given a vector.
- Handles the proportion-equals-0 and proportion-equals-1 boundary cases.
- Returns `NA`-filled output (with zero counts) for empty input data.
- Produces a wider interval at a higher confidence level.
- Matches against the `S_MUT` mutation profile when `mut = TRUE`.

### Layer 3 — Integration tests

Slower tests that exercise the nowcast model functions (`svymultinom` and
`se.multinom`) on `example_data.RDS`. These check **structural invariants**
rather than exact numeric outputs, so they remain valid even when the
underlying data is refreshed. They confirm:

- `svymultinom` returns the expected list structure and a fitted `multinom`
  object; when the Hessian is invertible, the design-adjusted covariance
  matrix is square, symmetric, and has a non-negative diagonal.
- `se.multinom` returns predicted proportions that sum to 1 and lie within
  [0, 1].
- Standard errors from `se.multinom` are non-negative and finite.
- Composite (aggregated) variant proportions equal the sum of their
  components.
- Predictions change continuously across nearby time points.

These tests are skipped automatically if `example_data.RDS` is unavailable.

## How to run the tests

From an R session, set the working directory so that the required files (listed
above) are reachable, then run the test file:

```r
testthat::test_file("test/test_variant_surveillance.R")
```

Alternatively, from the command line:

```sh
Rscript -e 'testthat::test_file("test/test_variant_surveillance.R")'
```

If your repository is set up as an R package or uses `testthat`'s standard
directory layout, you can instead run the whole suite with:

```r
testthat::test_dir("test")
```

## Interpreting the output

`testthat` prints a summary as it runs. A passing run reports each test context
with no failures. A typical successful run looks like:

```
== Testing test_variant_surveillance.R ===========================
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 35 ]
```

- **PASS** — the test's expectations were met.
- **FAIL** — an expectation was not met; the output shows the expected and
  actual values and the line number.
- **SKIP** — the test was skipped (most commonly the Layer 3 integration tests
  when `example_data.RDS` is not present). Skips are expected in that case and
  are not errors.

The Layer 1 and Layer 2 tests should complete in well under a second. The Layer
3 integration tests take longer because they fit the multinomial model on real
data; the prepared survey objects are cached, so the cost is paid only once per
run.

## Notes

- The test script sources `variant_surveillance_system.R` to obtain the
  unchanged helper functions. Sourcing wraps the call so that any error from
  the data-loading step does not prevent the already-defined functions from
  being used. If you prefer a cleaner setup, consider factoring the function
  definitions out of both scripts into a shared `helpers.R` and sourcing that
  instead.
- The integration tests use a minimal version of the data-preparation steps
  and a reduced variant list for speed. They are intended to catch structural
  regressions, not to reproduce the full analysis.
