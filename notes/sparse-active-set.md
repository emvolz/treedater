# A sparse tree‑Laplacian active‑set solver for treedater

## Summary

`treedater`'s inner node‑time optimisation was `O(n³)` per iteration, which
dominated run time on large trees. This change adds a **sparse Schur‑complement
active‑set** solver that is `O(n · #active‑constraints)` per solve — many times
faster on large trees — while returning numerically identical estimates. It is
exposed through a new value of the existing `clsSolver` argument,
`clsSolver = "sparse"`, which is now the **default**. The previous dense solver
remains available as `clsSolver = "limSolve"`.

| n (tips) | `limSolve` (dense) | `sparse` (new) | speed‑up | Δ tMRCA |
|---------:|-------------------:|---------------:|---------:|--------:|
|     250  |     0.083 s        |    0.019 s     |     4×   | 2.5e‑12 |
|     500  |     1.045 s        |    0.022 s     |    47×   | 1.4e‑12 |
|   1 000  |    14.3   s        |    0.067 s     |   214×   | 3.4e‑12 |
|   2 000  |    45.2   s        |    0.183 s     |   247×   | 3.6e‑12 |
|   4 000  |   (≈ 10 min)       |    0.516 s     |    —     |    —    |
|   8 000  |   (impractical)    |    1.395 s     |    —     |    —    |

(strict clock, single start condition; random binary trees with exponential
branch lengths and clock‑like tip dates. Δ tMRCA is the difference between the
sparse and dense estimates — i.e. they agree to `~1e‑12`.)

## The bottleneck

For a fixed set of per‑branch rate multipliers `omega`, treedater estimates the
`p = n − 1` internal‑node times by solving a weighted least‑squares problem
subject to node‑ordering (temporal) constraints — every parent must be at least
as old as each of its children/tips:

```
minimise    || sqrt(W) · (omega · A0 · x − B) ||²
subject to  Ain · x ≤ bin           (t_parent ≤ t_child ,  t_parent ≤ t_tip)
```

This quadratic program is re‑solved on **every** coordinate‑descent iteration of
the clock fit. The original solver, `.optim.Ti5.constrained.limsolve`, handed it
to `limSolve::lsei`, which forms the dense `p × p` normal matrix `AᵀWA` and runs
a dense active‑set QP — both steps are `O(n³)` (with an `O(n²)` memory
footprint). For `n` in the hundreds this already dominates; beyond ~2 000 tips it
is impractical.

## The key structure

Two observations make the problem sparse:

1. **`AᵀWA` is a weighted tree Laplacian.** Each row of `A0` has at most two
   non‑zeros (`−1` at the parent column, `+1` at the internal‑child column), so
   `AᵀWA` has the sparsity of the tree: a diagonal plus one off‑diagonal entry
   per internal edge (`≈ 2n` non‑zeros, not `n²`). Tip edges contribute only to
   the parent's diagonal, which "grounds" the Laplacian and makes it symmetric
   positive‑definite.

2. **The constraints are single edge rows.** Each temporal constraint involves
   at most two node times (`x_parent − x_child ≤ 0`, or `x_parent ≤ t_tip`), so
   the active constraints form a small, sparse sub‑matrix of the tree's incidence
   matrix.

## The algorithm (`.optim.Ti5.sparse.activeset`)

We solve the constrained QP with a primal **active‑set** method whose linear
algebra is driven by a **single sparse factorisation** of the tree Laplacian
`L = AᵀWA`:

1. **Factor once.** Build `L` as a sparse matrix and factor it with a sparse
   Cholesky (`Matrix::Cholesky`, i.e. CHOLMOD). On a tree, the elimination tree
   has no fill‑in, so the factor is `O(n)` to compute and to apply. (This is the
   linear‑algebra equivalent of a post‑order Gaussian elimination that walks the
   tree from tips to root.)
2. **Unconstrained solve.** `x0 = L⁻¹ c`, where `c = AᵀWB` is the right‑hand
   side (`O(n)`).
3. **Active‑set iterations.** Maintaining a working set of binding constraints
   `E x = f`:
   - apply `L⁻¹` to the active constraint columns, `Z = L⁻¹ Eᵀ` (each column
     `O(n)`; batched into one sparse solve);
   - form and solve the small dense Schur system `S μ = E x0 − f` with
     `S = E Z` (`|active| × |active|`) for the multipliers `μ`;
   - update `x = x0 − Z μ`;
   - **add** the most‑violated inactive constraint, or **release** the
     constraint with the most‑negative multiplier, and repeat.
4. **Terminate** when no constraint is violated and all multipliers are
   non‑negative (KKT satisfied).

The total cost is `O(n · #active)` — the tree factor is reused across all
active‑set iterations rather than refactoring a dense matrix each time. In
practice only a handful of temporal constraints bind, so the solver is close to
linear in `n`.

## The API option

`clsSolver` already selected the constrained‑least‑squares backend
(`"limSolve"` or `"mgcv"`); the new solver is a third choice and the default:

```r
clsSolver = c("sparse", "limSolve", "mgcv")     # "sparse" is the default
```

- `"sparse"` — the new tree‑Laplacian active‑set (requires the **Matrix**
  package, which is a "recommended" package shipped with R).
- `"limSolve"` — the previous dense `lsei` solver, unchanged (use this to
  reproduce pre‑change results exactly, or as a reference).
- `"mgcv"` — the `mgcv::pcls` path, unchanged.

The `"sparse"` path is defensive: if the **Matrix** package is unavailable, the
Cholesky factor is singular, or the active set fails to converge, it silently
falls back to `"limSolve"` (and then to `"mgcv"`), so results are always
produced.

Examples:

```r
dtr <- dater(tre, sts, s = 1500)                          # sparse (default), fast
dtr <- dater(tre, sts, s = 1500, clsSolver = "limSolve")  # dense reference
```

## Correctness / validation

- **Sparse vs. dense parity.** Across the phylodate test corpus and the bundled
  H3N2 example, `clsSolver = "sparse"` reproduces `clsSolver = "limSolve"`:
  - strict clock (deterministic): tMRCA and rate agree to `~1e‑12`
    (essentially bit‑for‑bit);
  - relaxed clocks (`uncorrelated`/`additive`, which wrap the solver in a
    Nelder‑Mead search): agree to `~1e‑5`, far tighter than any reporting
    precision — the tiny difference is the Nelder‑Mead trajectory amplifying a
    `~1e‑13` difference in the inner solve.
- **H3N2 example (177 tips).** Default and `"limSolve"` give identical
  `tMRCA = 1980.472`, `rate = 3.194e‑3`; downstream `outlierTips` and `parboot`
  run unchanged through the sparse path.
- **`R CMD check`** passes (examples OK, dependencies OK, Rd OK); the only NOTEs
  are pre‑existing/environmental (checking a git working tree; the existing
  vignette).

## Files changed

- `R/treedater0.R`
  - new `.optim.Ti5.sparse.activeset()` (the sparse Schur‑complement active‑set);
  - fit loop dispatches on `clsSolver == "sparse"` with graceful fallback;
  - `"sparse"` added as the default `clsSolver` choice in both `dater()` and the
    internal `.dater()`, and `clsSolver` is now forwarded from `dater()` to
    `.dater()` (it previously was not);
  - updated `@param clsSolver` documentation.
- `DESCRIPTION` — `Imports: Matrix`.
- `man/dater.Rd` — regenerated for the new `clsSolver` documentation.
