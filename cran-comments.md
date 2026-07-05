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
* win-builder, R-devel (2026-07-04 r90207 ucrt): **Status: OK**
* win-builder, R-release (R 4.6.1): **Status: OK**
* macOS builder: the service was unavailable at submission time (retry pending)

## R CMD check results

Both win-builder checks (R-devel and R-release) return **Status: OK** --
0 errors, 0 warnings, 0 notes.

The local `R CMD check --as-cran` is likewise clean apart from two artefacts of
optional external tools missing on the local machine, which do not occur on
CRAN / win-builder:

* WARNING: `'qpdf' is needed for checks on size reduction of PDFs` -- `qpdf`
  is not installed on the local machine.
* NOTE: `checking HTML version of manual ... Skipping checking HTML validation:
  no command 'tidy' found` -- HTML Tidy is not installed on the local machine.

The package's own examples, tests (`testthat`, edition 3) and vignette all build
and pass.
