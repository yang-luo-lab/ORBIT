# ====================================================================
# Tests for ORBIT_Rank
# ====================================================================
# ---------------------------------------------------------------------
# A. Basic structure of the return value
# ---------------------------------------------------------------------
test_that("ORBIT_Rank returns a data frame with the expected columns", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, rho = 0, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(ORBIT_Rank(sim$omics_list, direction, seed = 1))

  expect_s3_class(res, "data.frame")
  expect_named(res, c("Feature", "N", "P", "Direction"))
  expect_equal(nrow(res), 100)
  expect_setequal(res$Feature, sim$omics_list[[1]]$feature)
  expect_type(res$N, "integer")
  expect_type(res$P, "double")
  expect_type(res$Direction, "integer")
})

test_that("P, N, Direction have valid ranges and types", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = 0.2, seed = 1)
  )
  expect_true(all(res$P >= 0 & res$P <= 1, na.rm = TRUE))
  expect_true(all(res$N >= 0 & res$N <= 3))
  expect_true(all(res$Direction %in% c(-1L, 1L, NA_integer_)))
})

test_that("All-positive sign and direction yield Direction = +1", {
  set.seed(1)
  feat <- paste0("f", 1:50); ocol <- paste0("o", 1:3)
  stat <- matrix(runif(150), 50, 3)
  direction <- setNames(rep(1, 3), ocol)

  omics <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign = 1, stat = stat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  res <- suppressMessages(ORBIT_Rank(omics, direction, seed = 1))
  # All signed scores positive everywhere -> D > 0 -> Direction = 1.
  expect_true(all(res$Direction == 1L))
})

test_that("Flipping all signs flips Direction", {
  set.seed(1)
  feat <- paste0("f", 1:50); ocol <- paste0("o", 1:3)
  stat <- matrix(runif(150), 50, 3)
  direction <- setNames(rep(1, 3), ocol)

  pos <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign =  1, stat = stat[, j],
               stringsAsFactors = FALSE)
  }), ocol)
  neg <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign = -1, stat = stat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  r_pos <- suppressMessages(ORBIT_Rank(pos, direction, seed = 1))
  r_neg <- suppressMessages(ORBIT_Rank(neg, direction, seed = 1))

  m <- match(r_pos$Feature, r_neg$Feature)
  expect_true(all(r_pos$Direction == -r_neg$Direction[m]))
  expect_equal(r_pos$P, r_neg$P[m])  # P depends on |D|, not its sign
})


# ---------------------------------------------------------------------
# B. Effect of rho on P-values
# ---------------------------------------------------------------------
test_that("Larger positive rho gives larger P-values for most features", {
  sim <- ORBIT_simulate_null(n = 300, K = 4, rho = 0.5, seed = 1)
  direction <- setNames(rep(1, 4), names(sim$omics_list))

  res0 <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = 0,   seed = 1))
  res5 <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = 0.5, seed = 1))

  m <- match(res5$Feature, res0$Feature)
  ok <- res5$N == 4 & !is.na(res5$P) & !is.na(res0$P[m]) & res0$P[m] < 1
  # rho = 0.5 with K = 4 -> sigma_eff = 1 + 3 * 0.5 = 2.5; P should
  # systematically inflate for the vast majority of features.
  expect_gt(mean(res5$P[ok] > res0$P[m][ok]), 0.85)
  expect_gt(median(res5$P[ok] - res0$P[m][ok]), 0.01)
})


# ---------------------------------------------------------------------
# C. Statistical sanity
# ---------------------------------------------------------------------
test_that("ORBIT_Rank detects consistent multi-omic signal", {
  sim <- ORBIT_simulate_signal(
    n = 1000, K = 3,
    pi_null = 0.7, pi_up = 0.15, pi_down = 0.15,
    mu_signal = 3,
    rho_null = 0,
    rho_signal = 0.8,
    seed = 1
  )
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  res <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = 0, seed = 1)
  )
  fdr <- p.adjust(res$P, method = "BH")
  n_sig <- sum(fdr < 0.05, na.rm = TRUE)
  # 300 true signal features with rho_signal = 0.8 cross-omic
  # concordance should yield substantial detection at FDR 5%.
  expect_gt(n_sig, 30)
})

test_that("Null-only input gives uniformly-distributed P-values", {
  sim <- ORBIT_simulate_null(n = 1000, K = 3, rho = 0, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = 0, seed = 1)
  )
  pp <- res$P[!is.na(res$P)]
  expect_true(all(pp >= 0 & pp <= 1))
  expect_gt(mean(pp > 0.5), 0.45)
  expect_lt(mean(pp > 0.5), 0.55)
})


# ---------------------------------------------------------------------
# D. rho boundary handling (scalar mode)
# ---------------------------------------------------------------------

test_that("Negative rho within the valid range is accepted (no warning)", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)   # k_max = 3, lb = -0.5
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  expect_no_warning(
    res <- suppressMessages(
      ORBIT_Rank(sim$omics_list, direction, rho = -0.3, seed = 1)
    )
  )
  expect_s3_class(res, "data.frame")
  pp <- res$P[!is.na(res$P)]
  expect_true(all(pp >= 0 & pp <= 1))
})

test_that("Negative rho shrinks null variance -> smaller P than rho = 0", {
  sim <- ORBIT_simulate_null(n = 300, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  res_zero <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho =  0,    seed = 1))
  res_neg  <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = -0.3, seed = 1))

  m <- match(res_neg$Feature, res_zero$Feature)
  ok <- res_neg$N == 3 &
    !is.na(res_neg$P) & !is.na(res_zero$P[m]) &
    res_zero$P[m] < 1
  expect_true(any(ok))
  expect_true(all(res_neg$P[ok] <= res_zero$P[m][ok] + 1e-10))
  expect_gt(mean(res_neg$P[ok] < res_zero$P[m][ok]), 0.95)
})

test_that("rho at or below -1/(k_max - 1) is rejected with a clear error", {
  # K = 3 -> k_max = 3 -> lower bound = -0.5
  sim3 <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  d3   <- setNames(rep(1, 3), names(sim3$omics_list))

  expect_error(
    ORBIT_Rank(sim3$omics_list, d3, rho = -0.5),
    "lower bound"
  )
  expect_error(
    ORBIT_Rank(sim3$omics_list, d3, rho = -0.7),
    "lower bound"
  )
  sim2 <- ORBIT_simulate_null(n = 50, K = 2, seed = 1)
  d2   <- setNames(rep(1, 2), names(sim2$omics_list))
  expect_error(
    ORBIT_Rank(sim2$omics_list, d2, rho = -1),
    "lower bound"
  )
  expect_no_warning(
    res <- suppressMessages(
      ORBIT_Rank(sim2$omics_list, d2, rho = -0.9, seed = 1)
    )
  )
  expect_s3_class(res, "data.frame")
})

test_that("Lower bound is data-driven by k_max", {
  sim2 <- ORBIT_simulate_null(n = 50, K = 2, seed = 1)
  d2   <- setNames(rep(1, 2), names(sim2$omics_list))
  expect_no_warning(
    suppressMessages(
      ORBIT_Rank(sim2$omics_list, d2, rho = -0.6, seed = 1)
    )
  )
  sim3 <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  d3   <- setNames(rep(1, 3), names(sim3$omics_list))
  expect_error(
    ORBIT_Rank(sim3$omics_list, d3, rho = -0.6),
    "lower bound"
  )
})

test_that("rho = 1 is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = 1),
    "strictly less than 1"
  )
})

test_that("rho > 1 is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = 1.5),
    "must not exceed 1"
  )
})

test_that("Non-numeric or non-finite rho is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = NA),
    "finite number"
  )
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = "0.3"),
    "finite number"
  )
})


# ---------------------------------------------------------------------
# E. stat_direction
# ---------------------------------------------------------------------
test_that("stat_direction = 'larger_better' runs and produces valid P", {
  set.seed(7)
  feat <- paste0("feature", 1:60); ocol <- paste0("omic", 1:3)
  rank_stat <- abs(matrix(rnorm(60 * 3), 60, 3))
  sign_mat  <- matrix(sample(c(-1, 1), 60 * 3, replace = TRUE), 60, 3)
  direction <- setNames(rep(1, 3), ocol)

  omics <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat,
               sign    = sign_mat[, j],
               stat    = rank_stat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  res <- suppressMessages(
    ORBIT_Rank(omics, direction,
               stat_direction = "larger_better", seed = 1)
  )
  expect_true(all(res$P >= 0 & res$P <= 1, na.rm = TRUE))
})

test_that("stat_direction flips which end is treated as significant", {
  set.seed(42)
  feat <- paste0("f", 1:200); ocol <- paste0("o", 1:2)
  stat <- matrix(runif(400), 200, 2)
  direction <- setNames(c(1, 1), ocol)
  omics <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign = 1, stat = stat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  r_small <- suppressMessages(
    ORBIT_Rank(omics, direction,
               stat_direction = "smaller_better", seed = 1))
  r_large <- suppressMessages(
    ORBIT_Rank(omics, direction,
               stat_direction = "larger_better", seed = 1))

  ms <- match(feat, r_small$Feature)
  ml <- match(feat, r_large$Feature)
  expect_false(isTRUE(all.equal(r_small$P[ms], r_large$P[ml])))
})


# ---------------------------------------------------------------------
# F. Input validation (errors)
# ---------------------------------------------------------------------
test_that("Non-list / unnamed omics_list is rejected", {
  expect_error(
    ORBIT_Rank(matrix(1, 3, 3), direction = c(a = 1)),
    "named list"
  )
  unnamed <- list(
    data.frame(feature = "f1", sign = 1, stat = 0.5,
               stringsAsFactors = FALSE),
    data.frame(feature = "f1", sign = 1, stat = 0.5,
               stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_Rank(unnamed, direction = c(1, 1)),
    "named list"
  )
})

test_that("Non-data-frame omic element is rejected", {
  omics <- list(a = matrix(1, 3, 3))
  expect_error(
    ORBIT_Rank(omics, c(a = 1)),
    "not a data frame"
  )
})

test_that("Missing required column is rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2"),
                   stat = c(0.1, 0.2),  # no `sign`
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_Rank(omics, c(a = 1)),
    "missing column"
  )
})

test_that("Negative stat values are rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- sim$omics_list
  omics[[1]]$stat[1] <- -0.1
  direction <- setNames(rep(1, 3), names(omics))

  expect_error(
    ORBIT_Rank(omics, direction),
    "non-negative"
  )
})

test_that("NA in stat is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- sim$omics_list
  omics[[1]]$stat[1] <- NA
  direction <- setNames(rep(1, 3), names(omics))

  expect_error(
    ORBIT_Rank(omics, direction),
    "missing value"
  )
})

test_that("Bad sign values are rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 0, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_Rank(omics, c(a = 1)),
    "must be -1 or \\+1"
  )
})

test_that("NA in sign is rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, NA, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_Rank(omics, c(a = 1)),
    "missing value"
  )
})

test_that("Duplicate feature names within an omic are rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f1", "f2"),
                   sign    = c(1, 1, 1),
                   stat    = c(0.1, 0.2, 0.3),
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_Rank(omics, c(a = 1)),
    "duplicate"
  )
})

test_that("Unnamed direction is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  expect_error(
    ORBIT_Rank(sim$omics_list, c(1, 1, 1)),
    "named numeric"
  )
})

test_that("direction names must match omic names", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  expect_error(
    ORBIT_Rank(sim$omics_list, c(wrong = 1, names = 1, here = 1)),
    "must match names of"
  )
})

test_that("Non-+/-1 or NA direction values are rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  oc <- names(sim$omics_list)

  expect_error(
    ORBIT_Rank(sim$omics_list, setNames(c(1, 0, 1), oc)),
    "must contain only -1 or \\+1"
  )
  expect_error(
    ORBIT_Rank(sim$omics_list, setNames(c(1, 2, 1), oc)),
    "must contain only -1 or \\+1"
  )
  expect_error(
    ORBIT_Rank(sim$omics_list, setNames(c(1, NA, 1), oc)),
    "must contain only -1 or \\+1"
  )
})


# ---------------------------------------------------------------------
# G. Feature membership across omics
# ---------------------------------------------------------------------
test_that("Features missing from some omics get reduced N", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  omics <- sim$omics_list
  feat  <- omics[[1]]$feature

  omics[[1]] <- omics[[1]][!omics[[1]]$feature %in% feat[1:5], , drop = FALSE]
  omics[[2]] <- omics[[2]][!omics[[2]]$feature %in% feat[10:12], , drop = FALSE]

  direction <- setNames(rep(1, 3), names(omics))
  res <- suppressMessages(ORBIT_Rank(omics, direction, seed = 1))

  m <- match(feat, res$Feature)
  expect_equal(res$N[m[1:5]],   rep(2L, 5))
  expect_equal(res$N[m[10:12]], rep(2L, 3))
  expect_equal(res$N[m[20:30]], rep(3L, 11))
})


# ---------------------------------------------------------------------
# H. direction order-independence
# ---------------------------------------------------------------------
test_that("direction can be supplied in any order", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  oc  <- names(sim$omics_list)
  direction <- setNames(c(1, 1, -1), oc)

  res1 <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction,             seed = 42))
  res2 <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction[c(3, 1, 2)], seed = 42))
  expect_equal(res1, res2)
})


# ---------------------------------------------------------------------
# I. seed parameter
# ---------------------------------------------------------------------
test_that("Same seed yields reproducible output", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  out1 <- suppressMessages(ORBIT_Rank(sim$omics_list, direction, seed = 42))
  out2 <- suppressMessages(ORBIT_Rank(sim$omics_list, direction, seed = 42))
  expect_identical(out1, out2)
})

test_that("seed does not pollute the global RNG state", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  set.seed(123)
  before <- runif(1)

  set.seed(123)
  suppressMessages(ORBIT_Rank(sim$omics_list, direction, seed = 999))
  after <- runif(1)

  expect_identical(before, after)
})

test_that("Bad seed values are rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  expect_error(
    ORBIT_Rank(sim$omics_list, direction, seed = NA),
    "single finite number"
  )
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, seed = c(1, 2)),
    "single finite number"
  )
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, seed = "a"),
    "single finite number"
  )
})


# ---------------------------------------------------------------------
# J. Matrix rho mode
# ---------------------------------------------------------------------
# Helper: build a constant-off-diagonal rho matrix
make_const_rho_mat <- function(omic_names, val) {
  K <- length(omic_names)
  m <- matrix(val, K, K, dimnames = list(omic_names, omic_names))
  diag(m) <- 1
  m
}

test_that("Matrix rho with constant off-diagonal equals scalar rho", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  oc <- names(sim$omics_list)

  rho_val <- 0.3
  rho_mat <- make_const_rho_mat(oc, rho_val)

  res_scalar <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = rho_val, seed = 1))
  res_matrix <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = rho_mat, seed = 1))

  expect_equal(res_scalar$P, res_matrix$P, tolerance = 1e-10)
  expect_equal(res_scalar$Direction, res_matrix$Direction)
})

test_that("Matrix rho is invariant to row/col permutation", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  oc <- names(sim$omics_list)

  rho_mat <- matrix(c(1.0, 0.20, 0.10,
                      0.20, 1.0, 0.30,
                      0.10, 0.30, 1.0),
                    nrow = 3, dimnames = list(oc, oc))
  rho_perm <- rho_mat[c(3, 1, 2), c(3, 1, 2)]   # same matrix, reordered

  res1 <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = rho_mat,  seed = 1))
  res2 <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, rho = rho_perm, seed = 1))

  expect_equal(res1, res2)
})

test_that("Matrix rho gives per-feature eff_rho based on presence pattern", {
  feat <- paste0("f", 1:100); oc <- paste0("o", 1:3)
  set.seed(1)
  stat_mat <- matrix(runif(100 * 3), 100, 3)
  direction <- setNames(rep(1, 3), oc)

  omics <- setNames(list(
    data.frame(feature = feat, sign = 1, stat = stat_mat[, 1],
               stringsAsFactors = FALSE),
    data.frame(feature = feat, sign = 1, stat = stat_mat[, 2],
               stringsAsFactors = FALSE),
    # omic 3 only has the second half of features
    data.frame(feature = feat[51:100], sign = 1,
               stat = stat_mat[51:100, 3],
               stringsAsFactors = FALSE)
  ), oc)

  rho_mat <- matrix(c(1.0, 0.01, 0.90,
                      0.01, 1.0, 0.90,
                      0.90, 0.90, 1.0),
                    nrow = 3, dimnames = list(oc, oc))

  res <- suppressMessages(
    ORBIT_Rank(omics, direction, rho = rho_mat, seed = 1))

  m <- match(feat, res$Feature)
  expect_true(all(res$N[m[1:50]]   == 2L))
  expect_true(all(res$N[m[51:100]] == 3L))
  expect_true(all(res$P[m] >= 0 & res$P[m] <= 1, na.rm = TRUE))
})

test_that("Matrix rho with NA off-diagonal triggers warning and uses available pairs", {
  feat <- paste0("f", 1:100); oc <- paste0("o", 1:3)
  set.seed(1)
  stat_mat <- matrix(runif(100 * 3), 100, 3)
  direction <- setNames(rep(1, 3), oc)
  omics <- setNames(lapply(1:3, function(j)
    data.frame(feature = feat, sign = 1, stat = stat_mat[, j],
               stringsAsFactors = FALSE)
  ), oc)
  rho_mat <- matrix(c(1.0, 0.20, 0.10,
                      0.20, 1.0, NA,
                      0.10, NA,  1.0),
                    nrow = 3, dimnames = list(oc, oc))

  expect_warning(
    res <- suppressMessages(
      ORBIT_Rank(omics, direction, rho = rho_mat, seed = 1)),
    "NA pairs"
  )
  expect_true(all(res$P >= 0 & res$P <= 1, na.rm = TRUE))
})

test_that("Matrix rho rejects wrong dimensions", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  oc <- names(sim$omics_list)

  bad <- matrix(0.2, 2, 2, dimnames = list(oc[1:2], oc[1:2]))
  diag(bad) <- 1
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad),
    "omics"
  )

  bad2 <- matrix(0.2, 3, 4)
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad2),
    "square"
  )
})

test_that("Matrix rho rejects missing / wrong names", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  oc <- names(sim$omics_list)

  bad <- matrix(0.2, 3, 3); diag(bad) <- 1
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad),
    "rownames and colnames"
  )

  bad2 <- make_const_rho_mat(c("x", "y", "z"), 0.2)
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad2),
    "must match"
  )
})

test_that("Matrix rho rejects non-symmetric / bad diagonal / off-diag >= 1", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  oc <- names(sim$omics_list)

  bad1 <- matrix(c(1.0, 0.2, 0.1,
                   0.3, 1.0, 0.4,
                   0.1, 0.4, 1.0),
                 nrow = 3, dimnames = list(oc, oc))
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad1),
    "symmetric"
  )

  bad2 <- make_const_rho_mat(oc, 0.2)
  bad2[1, 1] <- 0.99
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad2),
    "diagonal"
  )

  bad3 <- make_const_rho_mat(oc, 1.0)
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad3),
    "strictly less than 1"
  )

  bad4 <- make_const_rho_mat(oc, 1.5)
  expect_error(
    ORBIT_Rank(sim$omics_list, direction, rho = bad4),
    "strictly less than 1"
  )
})

test_that("Matrix rho with non-positive effective inflation: P = NA + warning", {
  feat <- paste0("f", 1:50); oc <- paste0("o", 1:4)
  set.seed(1)
  stat_mat <- matrix(runif(50 * 4), 50, 4)
  direction <- setNames(rep(1, 4), oc)

  omics <- setNames(lapply(1:4, function(j)
    data.frame(feature = feat, sign = 1, stat = stat_mat[, j],
               stringsAsFactors = FALSE)
  ), oc)

  rho_mat <- make_const_rho_mat(oc, -0.4)

  expect_warning(
    res <- suppressMessages(
      ORBIT_Rank(omics, direction, rho = rho_mat, seed = 1)),
    "non-positive"
  )
  expect_true(all(is.na(res$P)))
  expect_true(all(is.na(res$Direction)))
})
