# treedater 2.0.0

## New features

* **Sparse node-time solver (`clsSolver = "sparse"`), now the default.**
  The temporally-constrained least-squares step that estimates internal-node
  times is now solved with a sparse tree-Laplacian **Goldfarb–Idnani dual
  active-set** instead of forming and solving a dense quadratic program. Each
  step applies `L⁻¹` to one constraint column (an `O(n)` sparse solve) and updates
  a small dense Schur factor. End-to-end this is roughly 6× faster than the dense
  solver on a real 1,000-tip relaxed-clock fit and ~15× at 2,000 tips, while
  returning numerically identical estimates (tMRCA/rate agree with the dense
  solver to ~1e-8 or better).

  Because it is a *dual* method, it also converges on **degenerate** constraint
  sets — many linearly-dependent binding constraints, common on large,
  non-clock-like trees — by exchanging constraints, where the earlier primal
  active-set had to fall back to the dense `quadprog` solver (e.g. ~14× faster
  than that fallback at 4,000 tips).

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
