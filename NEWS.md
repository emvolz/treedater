# treedater (development version)

## New features

* **Sparse node-time solver (`clsSolver = "sparse"`), now the default.**
  The temporally-constrained least-squares step that estimates internal-node
  times is now solved with a sparse tree-Laplacian Schur-complement active-set
  method instead of forming and solving a dense quadratic program. It is
  `O(n · #active-constraints)` per solve rather than `O(n³)`, giving large
  speed-ups on big trees (e.g. ~200× at 2,000 tips; trees of 8,000+ tips that
  were previously impractical now fit in a couple of seconds) while returning
  numerically identical estimates (tMRCA/rate agree with the old solver to
  ~1e-12 under a strict clock).

  The previous dense solver is still available as `clsSolver = "limSolve"`, and
  `"mgcv"` is unchanged. The sparse path requires the **Matrix** package (a
  recommended package shipped with R) and falls back to `"limSolve"`
  automatically if it is unavailable or fails.

  See `notes/sparse-active-set.md` for details.

## Bug fixes

* `dater()` now forwards its `clsSolver` argument to the internal fitting
  routine (previously the top-level `clsSolver` argument was ignored).
