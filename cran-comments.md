## Submission

This is a major release (version 2.0.0) of an existing CRAN package
(current CRAN version 1.0.2).

* The internal temporally-constrained least-squares step that estimates node
  times is now solved with a sparse tree-Laplacian (Goldfarb-Idnani dual
  active-set) method instead of the dense `limSolve::lsei`. Estimates are
  numerically identical to the previous release (times and rates agree to
  ~1e-8 or better on real data); the change gives substantial speed-ups on
  large trees.

* Dependency change (the reason for the major version bump): **`limSolve` has
  been removed from `Depends`.** The dense reference / fallback solver now uses
  **`quadprog`** (added to `Imports`, together with `Matrix`, a "recommended"
  package). The `clsSolver = "limSolve"` argument value is accordingly replaced
  by `clsSolver = "quadprog"`.

* This package has **no reverse dependencies on CRAN** (checked with
  `tools::package_dependencies(reverse = TRUE)` against the current CRAN
  package database), so the dependency and argument changes do not affect any
  other package.

## Test environments

* Local: Ubuntu 24.04, R 4.5.0 (x86_64-pc-linux-gnu)
* win-builder (R-release and R-devel) -- to be run before submission
* macOS builder -- to be run before submission

## R CMD check results

`R CMD check --as-cran` reports **0 errors and 0 warnings**, and no notes other
than the two below, which are artefacts of optional external tools missing on
the local build machine and do not occur on CRAN's check systems:

* WARNING: `'qpdf' is needed for checks on size reduction of PDFs` -- `qpdf`
  is not installed on the local machine.
* NOTE: `checking HTML version of manual ... Skipping checking HTML validation:
  no command 'tidy' found` -- HTML Tidy is not installed on the local machine.

The package's own examples, tests (`testthat`, edition 3) and vignette all build
and pass.
