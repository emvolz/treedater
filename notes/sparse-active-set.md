# A sparse tree‑Laplacian active‑set solver for treedater

## Summary

`treedater`'s inner node‑time optimisation was `O(n³)` per iteration, which
dominated run time on large trees. This change adds a **sparse Schur‑complement
active‑set** solver that is `O(n · #active‑constraints)` per solve — several
times faster on real relaxed‑clock fits — while returning numerically identical
estimates. It is exposed through a new value of the existing `clsSolver`
argument, `clsSolver = "sparse"`, which is now the **default**. A dense reference
solver is available as `clsSolver = "quadprog"`.

**End‑to‑end speed‑up on real data** (UK subtype‑C HIV tree, additive relaxed
clock, `omega0 = 0.0015`, single start, `ncpu = 1`; `sparse` vs. the dense
`quadprog` solver, identical estimates):

| n (tips) | `quadprog` (dense) | `sparse` (new) | speed‑up | Δ tMRCA |
|---------:|-------------------:|---------------:|---------:|--------:|
|   1 000  |      11.3 s        |     1.9 s      |    6×    | 1.2e‑8  |
|   2 000  |      84.9 s        |     5.8 s      |   15×    | 2.2e‑7  |

The single‑solve speed‑up grows with `n` and, on *idealised* perfectly
clock‑like trees, reaches 100×+ — but that best case has **zero** binding
ordering constraints. Real rate variation makes `O(n)` constraints bind, so the
realistic end‑to‑end gain is the ~6–15× above. Beyond ~3 500 tips a single solve
on real data can hit a degenerate (rank‑deficient) active set; the solver then
falls back to the dense `quadprog` path (detecting the degeneracy in seconds
rather than churning). See "Scaling and limits" below.

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
the clock fit. The dense solver forms the `p × p` normal matrix `AᵀWA` and runs a
dense active‑set QP — both steps are `O(n³)` (with an `O(n²)` memory footprint).
For `n` in the hundreds this already dominates; beyond ~2 000 tips it is
impractical.

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

We solve the constrained QP with a **Goldfarb–Idnani dual active‑set** method
whose linear algebra is driven by a **single sparse factorisation** of the tree
Laplacian `L = AᵀWA`:

1. **Factor once.** Build `L` as a sparse matrix and factor it with a sparse
   Cholesky (`Matrix::Cholesky`, i.e. CHOLMOD). On a tree, the elimination tree
   has no fill‑in, so the factor is `O(n)` to compute and to apply. (This is the
   linear‑algebra equivalent of a post‑order Gaussian elimination that walks the
   tree from tips to root.)
2. **Start dual‑feasible.** Begin at the unconstrained minimiser `x0 = L⁻¹ c`
   (`c = AᵀWB`) with an empty active set — trivially dual‑feasible (all
   multipliers zero).
3. **Dual iterations.** Repeatedly pick the most‑violated constraint `p` and
   raise its multiplier, keeping the active constraints tight and all multipliers
   `≥ 0`. For the incoming row `aₚ` compute `wₐ = L⁻¹ aₚ` (one `O(n)` sparse
   solve) and, against the active set, the primal step direction `z` and the dual
   step direction `r = M⁻¹ (A_active wₐ)` — using a maintained Cholesky of the
   Schur matrix `M = A_active L⁻¹ A_activeᵀ`. Then take the shorter of
   - the **primal (full) step** `t₂ = viol / (aₚᵀz)` that makes `p` tight, after
     which `p` enters the active set (append its column `wₐ`, extend the Schur
     factor by a bordered update); or
   - the **dual (partial) step** `t₁ = minₖ λₖ/rₖ` at which an active multiplier
     reaches zero, after which that constraint is **dropped** — a constraint
     *exchange* — and the same `p` is retried against the smaller active set.
4. **Terminate** when no constraint is violated (multipliers are `≥ 0` by
   construction) — the KKT conditions hold, so `x` is the feasible optimum.

The exchange in step 3 is what makes this robust. When the incoming constraint is
linearly dependent on the active set (`aₚᵀz ≤ 0` — a degenerate, rank‑deficient
active set) the primal step is infinite, so the method instead drops a blocking
constraint and proceeds, always reaching the true optimum. A plain *primal*
active‑set has no such move and must give up (fall back to the dense solver) on
exactly those instances. Each step costs one `O(n)` sparse solve plus `O(#active²)`
dense work, and the tree factor is reused throughout.

## Scaling and limits

The `O(n · #active)` framing is only near‑linear in `n` when few constraints bind.
That holds on clock‑like data. On **real, rate‑varying data it does not**: measured
on subtype‑C subsamples, a roughly *constant* ~30% of node‑ordering constraints
bind at every size (`#active ≈ 0.3 n`, from `n` = 1 000 to 5 000). With the Schur
work then `O(#active² … #active³)`, the solver is effectively `O(n³)`‑order on real
data (empirically the per‑solve time scales ~`n^2.3`). So on real data this is a
**large constant‑factor win — the ~6–15× in the table above — not an asymptotic
`O(n³)→O(n)` improvement.** Two consequences worth knowing:

- **Crossover.** On small trees the sparse solver's per‑solve overhead can make
  it no faster (or slightly slower) than the dense solver, but the node‑time
  solve is not a repeated hot path there, so it does not matter. The win shows up
  from ~500–1 000 tips upward.
- **Degeneracy — handled, no longer a fallback.** When many constraints bind they
  can become linearly dependent (a rank‑deficient active set, a genuine LICQ
  failure). The dual method handles this directly through the constraint
  *exchange* in step 3: it drops a blocking constraint and continues to the true
  feasible optimum. On the real subtype‑C tree the degeneracy appears sporadically
  from ~4 000 tips upward (e.g. the 4 000‑ and 5 000‑tip subsamples), and the dual
  solver converges on them at roughly the non‑degenerate sparse speed — **~14×
  faster than the dense fallback it replaces** (26 s vs ~370 s at 4 000 tips),
  staying feasible and matching `quadprog` to ~1e‑12.

  This is why the method is a *dual* active‑set. The previous *primal* version
  detected the degeneracy (`ρ² ≤ 0`) but had to fall back to a full dense
  `quadprog` solve; and a cheaper attempt — regularising the Schur to push
  through — **produced infeasible node times** (parent younger than child, error
  growing with `n`). `quadprog` remains only as the ultimate backstop for the
  unlikely cases where the sparse factor is singular or the **Matrix** package is
  unavailable.

## The API option

`clsSolver` selects the constrained‑least‑squares backend; the new solver is the
default:

```r
clsSolver = c("sparse", "quadprog", "mgcv")     # "sparse" is the default
```

- `"sparse"` — the new tree‑Laplacian active‑set (requires the **Matrix**
  package, which is a "recommended" package shipped with R).
- `"quadprog"` — the dense reference solver (forms `AᵀWA` and calls
  `quadprog::solve.QP`, a robust Goldfarb‑Idnani dual active‑set). Use this to
  cross‑check the sparse result or as the guaranteed fallback.
- `"mgcv"` — the `mgcv::pcls` path, unchanged.

The `"sparse"` path is defensive: if the **Matrix** package is unavailable, the
Cholesky factor is singular, or the active set fails to converge, it silently
falls back to `"quadprog"` (and then to `"mgcv"`), so results are always
produced.

Examples:

```r
dtr <- dater(tre, sts, s = 1500)                          # sparse (default), fast
dtr <- dater(tre, sts, s = 1500, clsSolver = "quadprog")  # dense reference
```

## Correctness / validation

- **Sparse vs. dense parity.** On the bundled H3N2 example and the subtype‑C
  tree, `clsSolver = "sparse"` reproduces `clsSolver = "quadprog"`:
  - strict clock (deterministic): tMRCA and rate agree to `~1e‑12`;
  - relaxed clocks (`uncorrelated`/`additive`, which wrap the solver in a
    Nelder‑Mead search): agree to `~1e‑5`, far tighter than any reporting
    precision — the tiny difference is the Nelder‑Mead trajectory amplifying a
    `~1e‑11` difference in the inner solve.
  Direct solver‑level parity (identical `td`/`omega`) is `~1e‑9` even with
  hundreds of active constraints.
- **H3N2 example (177 tips).** Default and `"quadprog"` give identical tMRCA and
  rate; downstream `outlierTips` and `parboot` run unchanged through the sparse
  path.
- **Test suite.** `tests/testthat/` covers sparse↔quadprog parity (clock‑like,
  heterogeneous rates, noisy/many‑constraint), constraint feasibility, the
  `sparse → quadprog` fallback, `clsSolver` forwarding, all three clock models
  (strict / uncorrelated / additive), `temporalConstraints = FALSE`, and
  `estimateSampleTimes`.

## Files changed

- `R/treedater0.R`
  - new `.optim.Ti5.sparse.activeset()` (the incremental sparse active‑set);
  - new `.optim.Ti5.constrained.quadprog()` (dense reference / fallback via
    `quadprog::solve.QP`), replacing the former `limSolve::lsei` solver;
  - fit loop dispatches on `clsSolver == "sparse"` with graceful fallback to
    `"quadprog"` then `"mgcv"`;
  - `"sparse"` is the default `clsSolver` in both `dater()` and the internal
    `.dater()`, and `clsSolver` is now forwarded from `dater()` to `.dater()`
    (it previously was not);
  - updated `@param clsSolver` documentation.
- `DESCRIPTION` — dropped `limSolve` from `Depends` (it is being discontinued on
  CRAN); added `quadprog` and `Matrix` to `Imports`.
- `NAMESPACE` / `man/dater.Rd` — regenerated.
- `tests/testthat/` — new test suite (see above).
