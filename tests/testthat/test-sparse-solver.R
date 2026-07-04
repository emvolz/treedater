# Tests for the sparse tree-Laplacian active-set node-time solver
# (clsSolver = "sparse") added in PR #27.
#
# The sparse solver must return the SAME constrained-least-squares solution as
# the dense reference solver (.optim.Ti5.constrained.quadprog, i.e. quadprog::solve.QP),
# while satisfying the node-ordering (temporal) constraints Ain %*% x <= bin.

# --- fixture -----------------------------------------------------------------
# A random binary tree whose tip "sample times" are the root-to-tip distances.
# These are clock-consistent, so the unconstrained weighted least-squares
# solution is close to feasible and the sparse active-set converges (rather than
# falling back to the dense solver) -- which is what we want to exercise here.
make_td <- function(n, seed, s = 1000, cc = 1) {
  set.seed(seed)
  tr  <- ape::rtree(n)
  D   <- ape::node.depth.edgelength(tr)          # depth of every node from the root
  sts <- stats::setNames(D[seq_len(n)], tr$tip.label)  # tips are nodes 1..n
  treedater:::.make.tree.data(tr, sts, s = s, cc = cc)
}

# =============================================================================
test_that("sparse solver reproduces the dense quadprog solution (omega = 1)", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  for (n in c(25L, 50L, 100L)) for (seed in 1:3) {
    td <- make_td(n, seed)
    om <- rep(1, nrow(td$tre$edge))
    xs <- treedater:::.optim.Ti5.sparse.activeset(om, td)
    xq <- treedater:::.optim.Ti5.constrained.quadprog(om, td)
    label <- sprintf("n=%d seed=%d", n, seed)
    # sparse path must actually run (not silently fall back) on clock-like data
    expect_false(is.null(xs), info = paste("sparse returned NULL for", label))
    expect_equal(xs, xq, tolerance = 1e-8, info = label)
  }
})

# =============================================================================
test_that("sparse solver matches quadprog with heterogeneous rate multipliers", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  # relaxed-clock inner solves call the solver with per-branch omega != 1
  for (seed in 1:3) {
    td <- make_td(80L, seed)
    set.seed(100 + seed)
    om <- runif(nrow(td$tre$edge), 0.5, 2)
    xs <- treedater:::.optim.Ti5.sparse.activeset(om, td)
    xq <- treedater:::.optim.Ti5.constrained.quadprog(om, td)
    expect_false(is.null(xs), info = sprintf("seed=%d", seed))
    expect_equal(xs, xq, tolerance = 1e-7, info = sprintf("seed=%d", seed))
  }
})

# =============================================================================
test_that("sparse solver matches quadprog under many binding constraints (noisy data)", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  # Poisson-noised branch lengths make MANY node-ordering constraints bind, so the
  # incremental active-set has to add a long sequence of constraints (maintaining Z
  # column-by-column + the bordered-Cholesky Schur factor). This exercises that
  # machinery far beyond the near-feasible clock-like fixtures above.
  set.seed(11)
  tr <- ape::rtree(120L, br = function(k) rexp(k, 10))
  D  <- ape::node.depth.edgelength(tr)[seq_len(120L)]
  tr$edge.length <- rpois(nrow(tr$edge), tr$edge.length * 995) / 995
  td <- treedater:::.make.tree.data(tr, stats::setNames(D, tr$tip.label), s = 995, cc = 10)
  om <- rep(0.0015, nrow(td$tre$edge))
  xs <- treedater:::.optim.Ti5.sparse.activeset(om, td)
  skip_if(is.null(xs), "sparse fell back on this instance")
  xq <- treedater:::.optim.Ti5.constrained.quadprog(om, td)
  expect_equal(xs, xq, tolerance = 1e-6)
  expect_lte(max(as.numeric(td$Ain %*% xs) - td$bin), 1e-6)   # feasible
})

# =============================================================================
test_that("sparse solution satisfies the node-ordering (temporal) constraints", {
  skip_if_not_installed("Matrix")
  # Every parent must be no younger than its children/tips: Ain %*% x <= bin.
  for (n in c(30L, 75L)) for (seed in 1:3) {
    td <- make_td(n, seed)
    om <- rep(1, nrow(td$tre$edge))
    xs <- treedater:::.optim.Ti5.sparse.activeset(om, td)
    expect_false(is.null(xs))
    viol <- max(as.numeric(td$Ain %*% xs) - td$bin)
    expect_lte(viol, 1e-6)
  }
})

# =============================================================================
test_that("sparse falls back gracefully and dater() still agrees with quadprog", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  # If the sparse solver ever returns NULL, dater() must fall back to the dense
  # quadprog solver and produce the correct answer. Force that by mocking the
  # sparse solver to always return NULL, then compare against a genuine quadprog fit.
  set.seed(1)
  tr  <- ape::rtree(40L)
  sts <- stats::setNames(ape::node.depth.edgelength(tr)[1:40], tr$tip.label)

  set.seed(2)
  fit_quadprog <- dater(tr, sts, s = 1000, clock = "strict",
                        clsSolver = "quadprog", ncpu = 1, numStartConditions = 0)

  testthat::local_mocked_bindings(
    .optim.Ti5.sparse.activeset = function(omegas, td) NULL,
    .package = "treedater"
  )
  set.seed(2)
  fit_fallback <- dater(tr, sts, s = 1000, clock = "strict",
                        clsSolver = "sparse", ncpu = 1, numStartConditions = 0)

  expect_equal(fit_fallback$timeOfMRCA, fit_quadprog$timeOfMRCA, tolerance = 1e-8)
  expect_equal(unname(fit_fallback$Ti), unname(fit_quadprog$Ti), tolerance = 1e-8)
})

# =============================================================================
test_that("dater() forwards clsSolver to the inner solver (regression for PR #27 bug fix)", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  # On master, dater() accepted clsSolver but never passed it to .dater(), so the
  # sparse path could never be reached from dater(). Prove the argument is now
  # forwarded by observing whether the sparse solver is actually invoked.
  set.seed(1)
  tr  <- ape::rtree(40L)
  sts <- stats::setNames(ape::node.depth.edgelength(tr)[1:40], tr$tip.label)

  called <- new.env(parent = emptyenv())
  called$sparse <- FALSE
  testthat::local_mocked_bindings(
    .optim.Ti5.sparse.activeset = function(omegas, td) { called$sparse <- TRUE; NULL },
    .package = "treedater"
  )

  called$sparse <- FALSE
  set.seed(2)
  dater(tr, sts, s = 1000, clock = "strict", clsSolver = "quadprog",
        ncpu = 1, numStartConditions = 0)
  expect_false(called$sparse)   # quadprog path must NOT touch the sparse solver

  called$sparse <- FALSE
  set.seed(2)
  dater(tr, sts, s = 1000, clock = "strict", clsSolver = "sparse",
        ncpu = 1, numStartConditions = 0)
  expect_true(called$sparse)    # sparse path MUST reach it -> clsSolver forwarded
})

# =============================================================================
test_that("dater() default (sparse) matches clsSolver='quadprog' on the H3N2 example", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  tre <- ape::read.tree(system.file("extdata", "flu_h3n2_final_small.treefile",
                                    package = "treedater"))
  sts <- sampleYearsFromLabels(tre$tip.label, delimiter = "_")
  seqlen <- 1698

  for (clk in c("strict", "uncorrelated")) {
    set.seed(1)
    fs <- dater(tre, sts, s = seqlen, clock = clk, clsSolver = "sparse",
                ncpu = 1, numStartConditions = 0)
    set.seed(1)
    fq <- dater(tre, sts, s = seqlen, clock = clk, clsSolver = "quadprog",
                ncpu = 1, numStartConditions = 0)
    expect_equal(fs$timeOfMRCA, fq$timeOfMRCA, tolerance = 1e-5, info = clk)
    expect_equal(fs$mean.rate,  fq$mean.rate,  tolerance = 1e-5, info = clk)
    expect_equal(unname(fs$Ti), unname(fq$Ti),  tolerance = 1e-5, info = clk)
  }
})

# =============================================================================
test_that("dater() sparse matches quadprog under the additive clock", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  tre <- ape::read.tree(system.file("extdata", "flu_h3n2_final_small.treefile",
                                    package = "treedater"))
  sts <- sampleYearsFromLabels(tre$tip.label, delimiter = "_")
  set.seed(1)
  fs <- dater(tre, sts, s = 1698, clock = "additive", omega0 = 0.003,
              clsSolver = "sparse",   ncpu = 1, numStartConditions = 0)
  set.seed(1)
  fq <- dater(tre, sts, s = 1698, clock = "additive", omega0 = 0.003,
              clsSolver = "quadprog", ncpu = 1, numStartConditions = 0)
  expect_equal(fs$timeOfMRCA, fq$timeOfMRCA, tolerance = 1e-5)
  expect_equal(fs$mean.rate,  fq$mean.rate,  tolerance = 1e-5)
  expect_equal(unname(fs$Ti), unname(fq$Ti),  tolerance = 1e-5)
})

# =============================================================================
test_that("dater() runs with temporalConstraints=FALSE (unconstrained path)", {
  skip_if_not_installed("Matrix")
  # temporalConstraints=FALSE bypasses the constrained solver entirely (it uses the
  # unconstrained .optim.Ti0/lm path), so it must (a) still run after limSolve removal
  # and (b) be unaffected by clsSolver.
  set.seed(1)
  tr  <- ape::rtree(40L)
  sts <- stats::setNames(ape::node.depth.edgelength(tr)[1:40], tr$tip.label)
  set.seed(2)
  f  <- dater(tr, sts, s = 1000, clock = "strict", temporalConstraints = FALSE,
              clsSolver = "sparse",   ncpu = 1, numStartConditions = 0)
  set.seed(2)
  fq <- dater(tr, sts, s = 1000, clock = "strict", temporalConstraints = FALSE,
              clsSolver = "quadprog", ncpu = 1, numStartConditions = 0)
  expect_s3_class(f, "treedater")
  expect_true(is.finite(f$timeOfMRCA))
  expect_equal(unname(f$Ti), unname(fq$Ti))   # solver unused -> identical
})

# =============================================================================
test_that("dater() sparse matches quadprog with estimateSampleTimes", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  # estimateSampleTimes rewrites td$sts2 (the tip constraint bounds) every iteration,
  # so it feeds the constrained solver -- verify sparse and quadprog agree.
  tre <- ape::read.tree(system.file("extdata", "flu_h3n2_final_small.treefile",
                                    package = "treedater"))
  sts <- sampleYearsFromLabels(tre$tip.label, delimiter = "_")
  unc <- names(sts)[1:6]
  est <- data.frame(lower = sts[unc] - 1, upper = sts[unc] + 1, row.names = unc)
  set.seed(1)
  fs <- dater(tre, sts, s = 1698, clock = "strict", estimateSampleTimes = est,
              clsSolver = "sparse",   ncpu = 1, numStartConditions = 0)
  set.seed(1)
  fq <- dater(tre, sts, s = 1698, clock = "strict", estimateSampleTimes = est,
              clsSolver = "quadprog", ncpu = 1, numStartConditions = 0)
  expect_equal(fs$timeOfMRCA, fq$timeOfMRCA, tolerance = 1e-5)
  expect_equal(unname(fs$Ti), unname(fq$Ti), tolerance = 1e-5)
  expect_equal(unname(fs$sts[unc]), unname(fq$sts[unc]), tolerance = 1e-5)  # estimates agree
  # estimation actually ran and stayed within the specified bounds
  expect_true(all(fs$sts[unc] >= sts[unc] - 1 - 1e-6 & fs$sts[unc] <= sts[unc] + 1 + 1e-6))
})

# =============================================================================
test_that("outlierTips gives the same result from a sparse fit and a quadprog fit", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  skip_if_not_installed("quadprog")
  tre <- ape::read.tree(system.file("extdata", "flu_h3n2_final_small.treefile",
                                    package = "treedater"))
  sts <- sampleYearsFromLabels(tre$tip.label, delimiter = "_")
  set.seed(1)
  fs <- dater(tre, sts, s = 1698, clock = "strict", clsSolver = "sparse",
              ncpu = 1, numStartConditions = 0)
  set.seed(1)
  fq <- dater(tre, sts, s = 1698, clock = "strict", clsSolver = "quadprog",
              ncpu = 1, numStartConditions = 0)
  os <- outlierTips(fs)
  oq <- outlierTips(fq)
  expect_s3_class(os, "data.frame")
  expect_true(all(c("taxon", "q", "p") %in% colnames(os)))
  os <- os[order(os$taxon), ]; oq <- oq[order(oq$taxon), ]
  expect_equal(os$q, oq$q, tolerance = 1e-6)
  expect_equal(sum(os$q < 0.05), sum(oq$q < 0.05))   # same lineages flagged
})

# =============================================================================
test_that("parboot runs end-to-end through the sparse solver (smoke, tiny nreps)", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  # parboot re-fits nreps times, so keep nreps and the tree small: this is a
  # pipeline smoke test (the bootstrap refits go through the default sparse solver),
  # not a statistical check.
  set.seed(3)
  tr  <- ape::rtree(45L)
  sts <- stats::setNames(ape::node.depth.edgelength(tr)[1:45], tr$tip.label)
  set.seed(3)
  f <- dater(tr, sts, s = 1000, clock = "strict", clsSolver = "sparse",
             ncpu = 1, numStartConditions = 0)
  set.seed(3)
  pb <- parboot(f, nreps = 3, ncpu = 1, quiet = TRUE)
  expect_s3_class(pb, "bootTreedater")
  expect_length(pb$timeOfMRCA_CI, 2)
  expect_true(all(is.finite(pb$timeOfMRCA_CI)))
  expect_false(is.null(pb$meanRate_CI))
})

# =============================================================================
test_that(".optim.Ti0.sparse (unconstrained path) matches the dense lm solution", {
  skip_if_not_installed("Matrix")
  # the temporalConstraints=FALSE path uses .optim.Ti0, which now solves the
  # unconstrained WLS via the sparse tree Laplacian instead of a dense lm; check
  # the two agree (for both weighting modes).
  lm_ti0 <- function(omegas, td, svbr = FALSE) {
    A <- omegas * td$A0
    B <- td$B0
    B[td$tipEdges] <- td$B0[td$tipEdges] - unname(omegas[td$tipEdges] * td$sts2)
    w <- if (svbr) td$W0 / omegas else td$W0
    unname(stats::coef(stats::lm(B ~ A - 1, weights = w)))
  }
  for (seed in 1:3) for (svbr in c(FALSE, TRUE)) {
    td <- make_td(80L, seed)
    set.seed(seed + 50)
    om <- runif(nrow(td$tre$edge), 0.5, 2)
    xs <- treedater:::.optim.Ti0.sparse(om, td, svbr)
    xl <- lm_ti0(om, td, svbr)
    expect_false(is.null(xs), info = sprintf("seed=%d svbr=%s", seed, svbr))
    expect_equal(xs, xl, tolerance = 1e-6, info = sprintf("seed=%d svbr=%s", seed, svbr))
  }
})
