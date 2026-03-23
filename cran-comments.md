## Resubmission

This is a resubmission of v1.3.1, which was rejected due to excessive test
CPU time (22x CPU/elapsed ratio on Debian). Changes in v1.3.2:

* Limit RcppParallel threads to 2 during tests (CRAN policy compliance)
* Share MCMC fit fixtures across plot/table tests (reduces total fits from ~32 to ~12)
* Reduce MCMC iterations in structural tests from 3000 to 500
* Mark expensive recovery tests with skip_on_cran()
* Use configure script for platform-specific LAPACK linking (macOS vs Linux)

Test time is now ~22s CPU / ~14s elapsed (ratio 1.6x).

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.

* The DOI https://doi.org/10.1093/infdis/jiu186 (in CITATION and vignette)
  returns HTTP 403 to automated checks but resolves correctly in a browser.
  This is standard behaviour from Oxford Academic's bot detection.

## Test environments

* macOS (aarch64-apple-darwin20), R 4.2.0
* GitHub Actions: macOS-latest (R release), ubuntu-latest (R release, R devel), windows-latest (R release)
