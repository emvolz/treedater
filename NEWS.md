# treedater 2.0.0

## New features

* **Sparse node-time solver (`clsSolver = "sparse"`), now the default.**
  The temporally-constrained least-squares step that estimates internal-node
  times is now solved with a sparse tree-Laplacian Schur-complement active-set
  method instead of forming and solving a dense quadratic program. The active set
  is maintained incrementally (one sparse solve per binding constraint plus a
  bordered-Cholesky update of the Schur factor), giving `O(n · #active)` per solve
  rather than `O(n³)`. End-to-end this is roughly 6× faster than the dense solver
  on a real 1,000-tip relaxed-clock fit and ~15× at 2,000 tips, while returning
  numerically identical estimates (tMRCA/rate agree with the dense solver to
  ~1e-8 or better).

  The dense solver is available as `clsSolver = "quadprog"`, and `"mgcv"` is
  unchanged. The sparse path requires the **Matrix** package (a recommended
  package shipped with R) and falls back to `"quadprog"` automatically if it is
  unavailable or fails, so a result is always produced.

  See `notes/sparse-active-set.md` for details.

## Dependencies

* **`limSolve` is no longer required** (it is being discontinued on CRAN). The
  dense constrained-least-squares solver now uses **`quadprog`** (moved to
  `Imports`); `limSolve` has been dropped from `Depends`.

## Bug fixes

* `dater()` now forwards its `clsSolver` argument to the internal fitting
  routine (previously the top-level `clsSolver` argument was ignored).
