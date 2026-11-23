# Regression tests for the two "extra" mechanisms of the genome engine:
#   (1) non-allelic gene conversion, within- and between-chromosome
#       (p_conversion, conv_dist, conv_tract, p_conv_interchrom)
#   (2) local-duplication hotspot feedback (hotspot_alpha)
# Small/fast sims; each test targets one designed property of the model.

K       <- 3L
TARGETS <- rep(300L, 3L)
NGEN    <- 400L
MU      <- 1e-4
L       <- 178L

# number of distinct sequences carried by more than one chromosome
n_cross_chrom <- function(g)
  nrow(g[, .(nc = data.table::uniqueN(chrom)), by = bponly][nc > 1])

# size of the largest set of identical monomers within a single chromosome
max_family <- function(g) max(g[, .N, by = .(chrom, bponly)]$N)

# A genome with no duplication and no deletion: the governor's balance term is
# zero and each chromosome sits EXACTLY on its target, so any length change
# would have to come from conversion itself.
static_sim <- function(seed, ...) {
  set.seed(seed)
  sim_genome(K = K, target_sizes = TARGETS, n_generations = NGEN,
             mu_total = MU, init_l = L,
             p_local_dup = 0, p_distal_dup = 0, p_exchange = 0, ...)
}

# ---- gene conversion: the length contract --------------------------------

test_that("gene conversion rewrites sequences without changing monomer count", {
  plain <- static_sim(11, p_conversion = 0)
  conv  <- static_sim(11, p_conversion = 5e-3)

  # no duplication/deletion => both genomes sit exactly on target
  expect_identical(attr(plain, "array_sizes"), TARGETS)
  expect_identical(attr(conv,  "array_sizes"), TARGETS)
  expect_identical(nrow(conv), nrow(plain))

  # ... but the sequences are not the same ones
  expect_false(identical(conv$bponly, plain$bponly))

  # conversion homogenizes: fewer distinct monomer sequences survive
  expect_lt(data.table::uniqueN(conv$bponly),
            data.table::uniqueN(plain$bponly))
})

test_that("more conversion means more homogenization", {
  u <- vapply(c(0, 1e-3, 5e-3), function(p)
    data.table::uniqueN(static_sim(12, p_conversion = p)$bponly), 0L)
  expect_true(all(diff(u) < 0))
})

test_that("a longer conversion tract converts more monomers per event", {
  short <- static_sim(13, p_conversion = 2e-3,
                      conv_tract = list(type = "fixed", value = 1))
  long  <- static_sim(13, p_conversion = 2e-3,
                      conv_tract = list(type = "fixed", value = 10))
  expect_identical(nrow(long), nrow(short))          # still no length change
  expect_lt(data.table::uniqueN(long$bponly),
            data.table::uniqueN(short$bponly))
})

# ---- inter-chromosomal conversion ----------------------------------------

test_that("p_conv_interchrom > 0 creates cross-chromosome identical monomers", {
  # p_exchange = 0, so conversion is the ONLY route to cross-chrom sharing
  within <- static_sim(14, p_conversion = 5e-3, p_conv_interchrom = 0)
  across <- static_sim(14, p_conversion = 5e-3, p_conv_interchrom = 1)

  expect_identical(attr(across, "n_exchange"), 0)   # no block was translocated
  expect_identical(n_cross_chrom(within), 0L)
  expect_gt(n_cross_chrom(across), 0L)
  expect_identical(nrow(across), nrow(within))      # still no length change
})

test_that("inter-chromosomal conversion scatters copies rather than moving blocks", {
  # conversion transfers short tracts into occupied positions, so its shared
  # monomers do not come in the contiguous runs an exchanged block leaves behind
  mean_shared_run <- function(g) {
    shared <- g[, .(nc = data.table::uniqueN(chrom)), by = bponly][nc > 1, bponly]
    runs <- g[, {
      r <- rle(bponly %in% shared)
      .(len = r$lengths[r$values])
    }, by = chrom]$len
    mean(runs)
  }
  by_conv <- lapply(1:3, function(s) {
    set.seed(s)
    sim_genome(K = K, target_sizes = rep(400L, K), n_generations = NGEN,
               mu_total = MU, init_l = L, p_exchange = 0,
               p_conversion = 5e-3, p_conv_interchrom = 1)
  })
  by_xch <- lapply(1:3, function(s) {
    set.seed(s)
    sim_genome(K = K, target_sizes = rep(400L, K), n_generations = NGEN,
               mu_total = MU, init_l = L,
               p_exchange = 0.2, exchange_chunk_mean = 200)
  })
  expect_lt(mean(vapply(by_conv, mean_shared_run, 0.0)),
            mean(vapply(by_xch,  mean_shared_run, 0.0)))
})

# ---- hotspot feedback ------------------------------------------------------

# Regime chosen so identical families are big enough to see: fast local
# duplication of single monomers, mutation slow enough that families persist.
HOT_SEEDS <- 1:5

hot_sim <- function(alpha, seed) {
  set.seed(seed)
  sim_genome(K = 2L, target_sizes = rep(400L, 2L), n_generations = 400L,
             mu_total = 1e-4, init_l = L,
             p_local_dup = 5e-3,
             local_dist = list(type = "fixed", value = 1),
             del_dist   = list(type = "fixed", value = 1),
             hotspot_alpha = alpha)
}

hot_grid <- data.table::rbindlist(lapply(c(0, 3), function(a)
  data.table::rbindlist(lapply(HOT_SEEDS, function(s) {
    g <- hot_sim(a, s)
    data.table::data.table(alpha = a, seed = s,
                           max_fam = max_family(g), total = nrow(g))
  }))))

test_that("hotspot feedback builds larger identical families", {
  m0 <- mean(hot_grid[alpha == 0, max_fam])
  m3 <- mean(hot_grid[alpha == 3, max_fam])
  expect_gt(m3, 1.4 * m0)   # observed ratio ~2.2; margin is deliberately loose
})

test_that("hotspot feedback is mean-preserving in duplication mass", {
  # alpha rescales the per-monomer trigger probability so the expected number
  # of duplication events per chromosome is unchanged. The governor's deletion
  # rate is derived analytically from p_local_dup (not from realized events),
  # so the standing monomer count is a direct readout of realized duplication
  # mass: it must not shift when alpha turns on.
  t0 <- mean(hot_grid[alpha == 0, total])
  t3 <- mean(hot_grid[alpha == 3, total])
  expect_lt(abs(t3 - t0) / t0, 0.15)   # stochastic; observed ~0.5%
})

test_that("hotspot_alpha = 0 is the default and is validated", {
  expect_identical(eval(formals(sim_genome)$hotspot_alpha), 0)
  expect_identical(eval(formals(sim_genome)$p_conversion), 0)
  expect_identical(eval(formals(sim_genome)$p_conv_interchrom), 0)
  expect_error(sim_genome(K = 2L, target_sizes = c(50L, 50L),
                          n_generations = 10L, hotspot_alpha = -1),
               "hotspot_alpha")
})

# ---- reproducibility with both mechanisms on -------------------------------

test_that("same seed reproduces the identical genome with both mechanisms on", {
  args <- list(K = K, target_sizes = TARGETS, n_generations = NGEN,
               mu_total = MU, init_l = L,
               p_conversion = 3e-3, p_conv_interchrom = 0.5,
               conv_dist  = list(type = "gamma", shape = 1, scale = 20),
               conv_tract = list(type = "gamma", shape = 1, scale = 1),
               hotspot_alpha = 2)
  set.seed(99); a <- do.call(sim_genome, args)
  set.seed(99); b <- do.call(sim_genome, args)
  expect_identical(a, b)
})
