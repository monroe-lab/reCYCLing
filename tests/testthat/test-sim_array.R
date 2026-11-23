# Tests for the single-array simulator and the shared monomer data contract.

test_that("sim_array returns the documented structure", {
  set.seed(1)
  s <- sim_array(max_units = 300, init_l = 60, mu_total = 1e-4)
  expect_named(s, c("monomers", "pairs", "L_vec", "H_vec", "N_vec"))
  expect_named(s$monomers, c("num", "bponly", "dir"))
  expect_true(nrow(s$monomers) >= 300)
  expect_true(all(s$monomers$dir %in% c("+", "-")))
  expect_true(all(c("bponly", "chrom1", "num1", "dir1",
                    "chrom2", "num2", "dir2") %in% names(s$pairs)))
})

test_that("compute_pairs = FALSE skips the pair table", {
  set.seed(2)
  s <- sim_array(max_units = 200, init_l = 60, mu_total = 1e-4,
                 compute_pairs = FALSE)
  expect_identical(nrow(s$pairs), 0L)
})

test_that("same seed reproduces the identical array (both backends)", {
  set.seed(3); a <- sim_array(max_units = 200, init_l = 60, mu_total = 1e-4)
  set.seed(3); b <- sim_array(max_units = 200, init_l = 60, mu_total = 1e-4)
  expect_identical(a, b)
  set.seed(4); r1 <- sim_array(max_units = 100, init_l = 40, mu_total = 1e-4,
                               backend = "r", max_generations = 5000)
  set.seed(4); r2 <- sim_array(max_units = 100, init_l = 40, mu_total = 1e-4,
                               backend = "r", max_generations = 5000)
  expect_identical(r1, r2)
})

test_that("optional mechanisms are off by default", {
  frm <- formals(sim_array)
  expect_identical(eval(frm$p_distal_dup), 0)
  expect_identical(eval(frm$p_del_chunk), 0)
  expect_identical(eval(frm$p_terminal_del), 0)
  expect_identical(eval(frm$p_conversion), 0)
})

test_that("gene conversion rewrites monomers without changing array length", {
  # no duplication and no deletion: only conversion (and mutation) can act, so
  # the array length is fixed by construction and any change would be a bug
  static <- function(p_conversion) {
    set.seed(6)
    sim_array(init_l = 60, init_k0 = 40, max_generations = 600,
              p_local_dup = 0, p_del_chunk = 0, mu_total = 3e-4,
              p_conversion = p_conversion, compute_pairs = FALSE)
  }
  plain <- static(0)
  conv  <- static(0.01)
  expect_identical(nrow(plain$monomers), 40L)
  expect_identical(nrow(conv$monomers), nrow(plain$monomers))
  expect_false(identical(conv$monomers$bponly, plain$monomers$bponly))
  # conversion homogenizes: fewer distinct sequences at the same array length
  expect_lt(data.table::uniqueN(conv$monomers$bponly),
            data.table::uniqueN(plain$monomers$bponly))
})

test_that("gene conversion requires the C++ backend", {
  expect_error(
    sim_array(max_units = 100, init_l = 40, p_conversion = 0.01, backend = "r"),
    "C\\+\\+ backend")
})

test_that("pairs_identical works on the genome contract and finds cross-chrom pairs", {
  set.seed(5)
  g <- sim_genome(K = 3L, target_sizes = rep(200L, 3L), n_generations = 300L,
                  mu_total = 1e-4, init_l = 60L,
                  p_exchange = 0.3, exchange_chunk_mean = 50)
  ps <- pairs_identical(g)
  expect_true(all(c("chrom1", "chrom2") %in% names(ps)))
  expect_gt(nrow(ps[chrom1 != chrom2]), 0L)
})

test_that("pairs_identical does not modify its input", {
  dt <- data.table::data.table(bponly = c("AA", "AA", "CC"),
                               num = 1:3, dir = "+")
  before <- data.table::copy(dt)
  invisible(pairs_identical(dt))
  expect_identical(dt, before)
})
