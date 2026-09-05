# ====================================================================
# Tests for ORBIT_P (updated for current API)
# ====================================================================


# ---------------------------------------------------------------------
# A. Basic structure of the return value
# ---------------------------------------------------------------------
test_that("ORBIT_P returns a data frame with the expected columns", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, rho = 0, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(ORBIT_P(sim$omics_list, direction, rho = 0))

  expect_s3_class(res, "data.frame")
  expect_named(res, c("Feature", "N", "P", "Direction"))
  expect_equal(nrow(res), 100)
  expect_setequal(res$Feature, sim$omics_list[[1]]$feature)
  expect_type(res$N, "integer")
  expect_type(res$P, "double")
  expect_type(res$Direction, "integer")
})

test_that("Return value has a `rho` attribute echoing the rho used", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  # Scalar supplied -> scalar attribute
  res_s <- suppressMessages(
    ORBIT_P(sim$omics_list, direction, rho = 0.2))
  expect_equal(attr(res_s, "rho"), 0.2)

  # Auto-estimated -> matrix attribute
  res_a <- suppressMessages(ORBIT_P(sim$omics_list, direction))
  rho_a <- attr(res_a, "rho")
  expect_true(is.matrix(rho_a))
  expect_equal(dim(rho_a), c(3, 3))
  expect_true(all(abs(diag(rho_a) - 1) < 1e-8))
  expect_true(isSymmetric(rho_a))
})

test_that("P, N, Direction have valid ranges and types", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(
    ORBIT_P(sim$omics_list, direction, rho = 0.2)
  )
  expect_true(all(res$P >= 0 & res$P <= 1, na.rm = TRUE))
  expect_true(all(res$N >= 0 & res$N <= 3))
  expect_true(all(res$Direction %in% c(-1L, 1L, NA_integer_)))
})

test_that("All-positive sign and direction yield Direction = +1", {
  set.seed(1)
  feat <- paste0("f", 1:50); ocol <- paste0("o", 1:3)
  pmat <- matrix(runif(150), 50, 3)
  direction <- setNames(rep(1, 3), ocol)

  omics <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign = 1, stat = pmat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  res <- suppressMessages(ORBIT_P(omics, direction, rho = 0))
  expect_true(all(res$Direction == 1L))
})

test_that("Flipping all signs flips Direction; P unchanged", {
  set.seed(1)
  feat <- paste0("f", 1:50); ocol <- paste0("o", 1:3)
  pmat <- matrix(runif(150), 50, 3)
  direction <- setNames(rep(1, 3), ocol)

  pos <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign =  1, stat = pmat[, j],
               stringsAsFactors = FALSE)
  }), ocol)
  neg <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign = -1, stat = pmat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  r_pos <- suppressMessages(ORBIT_P(pos, direction, rho = 0))
  r_neg <- suppressMessages(ORBIT_P(neg, direction, rho = 0))

  m <- match(r_pos$Feature, r_neg$Feature)
  expect_true(all(r_pos$Direction == -r_neg$Direction[m]))
  expect_equal(r_pos$P, r_neg$P[m])
})


# ---------------------------------------------------------------------
# B. Auto-rho estimation (matrix mode)
# ---------------------------------------------------------------------
test_that("rho is auto-estimated as a matrix when not supplied", {
  sim <- ORBIT_simulate_null(n = 500, K = 3, rho = 0.4, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  expect_message(
    ORBIT_P(sim$omics_list, direction),
    "Auto-estimated rho matrix"
  )

  res <- suppressMessages(ORBIT_P(sim$omics_list, direction))
  rho_used <- attr(res, "rho")
  expect_true(is.matrix(rho_used))
})

test_that("Supplied rho suppresses the auto-estimation message", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.4, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  msgs <- capture_messages(
    ORBIT_P(sim$omics_list, direction, rho = 0.2)
  )
  expect_false(any(grepl("Auto-estimated", msgs)))
})

test_that("Independent null: auto-estimated rho matrix has median near 0", {
  sim <- ORBIT_simulate_null(n = 3000, K = 3, rho = 0, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(ORBIT_P(sim$omics_list, direction))
  rho_used <- attr(res, "rho")
  expect_true(is.matrix(rho_used))

  ut_vals <- rho_used[upper.tri(rho_used)]
  ut_vals <- ut_vals[is.finite(ut_vals)]
  expect_true(length(ut_vals) > 0)
  expect_lt(abs(median(ut_vals)), 0.02)
})

test_that("Correlated null: auto-estimated rho matrix recovers the truth", {
  true_rho <- 0.4
  sim <- ORBIT_simulate_null(n = 3000, K = 3, rho = true_rho, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(ORBIT_P(sim$omics_list, direction))
  rho_used <- attr(res, "rho")
  expect_true(is.matrix(rho_used))

  ut_vals <- rho_used[upper.tri(rho_used)]
  ut_vals <- ut_vals[is.finite(ut_vals)]
  # With per-pair rerank on full overlap (no signal filtering),
  # estimate should be close to truth under pure null.
  expect_lt(abs(median(ut_vals) - true_rho), 0.03)
})

test_that("Single-omic input sets rho = 0 with a message", {
  set.seed(1)
  omics <- list(
    a = data.frame(feature = paste0("f", 1:20),
                   sign    = sample(c(-1, 1), 20, replace = TRUE),
                   stat    = runif(20),
                   stringsAsFactors = FALSE)
  )
  direction <- c(a = 1)
  expect_message(
    res <- ORBIT_P(omics, direction),
    "Single-omic"
  )
  expect_equal(attr(res, "rho"), 0)
  expect_true(all(res$P >= 0 & res$P <= 1, na.rm = TRUE))
})


# ---------------------------------------------------------------------
# C. Effect of rho on P-values
# ---------------------------------------------------------------------
test_that("Larger rho gives larger P-values for most features", {
  sim <- ORBIT_simulate_null(n = 300, K = 4, rho = 0.5, seed = 1)
  direction <- setNames(rep(1, 4), names(sim$omics_list))

  res0 <- suppressMessages(ORBIT_P(sim$omics_list, direction, rho = 0))
  res5 <- suppressMessages(ORBIT_P(sim$omics_list, direction, rho = 0.5))

  m <- match(res5$Feature, res0$Feature)
  ok <- res5$N == 4 & !is.na(res5$P) & !is.na(res0$P[m]) & res0$P[m] < 1
  # rho = 0.5 with K = 4 -> sigma_eff = 1 + 3 * 0.5 = 2.5; P should
  # systematically inflate for the vast majority of features.
  expect_gt(mean(res5$P[ok] > res0$P[m][ok]), 0.85)
  expect_gt(median(res5$P[ok] - res0$P[m][ok]), 0.01)
})


# ---------------------------------------------------------------------
# D. Statistical sanity
# ---------------------------------------------------------------------
test_that("ORBIT_P detects consistent multi-omic signal", {
  sim <- ORBIT_simulate_signal(
    n = 1000, K = 3,
    pi_null = 0.8, pi_up = 0.1, pi_down = 0.1,
    mu_signal = 3,
    rho_null = 0,
    rho_signal = 0.8,
    seed = 1
  )
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(ORBIT_P(sim$omics_list, direction))
  fdr <- p.adjust(res$P, method = "BH")
  n_sig <- sum(fdr < 0.05, na.rm = TRUE)
  # 200 true signal features with rho_signal = 0.8 cross-omic
  # concordance should give substantial detection at FDR 5%.
  expect_gt(n_sig, 50)
})

test_that("Null-only input gives uniformly-distributed P-values", {
  sim <- ORBIT_simulate_null(n = 1000, K = 3, rho = 0, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  res <- suppressMessages(ORBIT_P(sim$omics_list, direction, rho = 0))
  pp <- res$P[!is.na(res$P)]
  expect_true(all(pp >= 0 & pp <= 1))
  # Under correct calibration P ~ Uniform(0, 1); with n = 1000,
  # SE ~ 0.016, so fraction above 0.5 should be in [0.45, 0.55].
  expect_gt(mean(pp > 0.5), 0.45)
  expect_lt(mean(pp > 0.5), 0.55)
})


# ---------------------------------------------------------------------
# E. P-value range validation
# ---------------------------------------------------------------------
test_that("Negative P-values are rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- sim$omics_list
  omics[[1]]$stat[1] <- -0.01
  direction <- setNames(rep(1, 3), names(omics))

  expect_error(
    ORBIT_P(omics, direction),
    "in \\[0, 1\\]"
  )
})

test_that("P > 1 is rejected with a hint to use ORBIT_Rank", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- sim$omics_list
  omics[[1]]$stat[1] <- 1.5
  direction <- setNames(rep(1, 3), names(omics))

  expect_error(
    ORBIT_P(omics, direction),
    "ORBIT_Rank"
  )
})

test_that("P = 0 is clamped with a warning", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- sim$omics_list
  omics[[1]]$stat[1] <- 0
  omics[[2]]$stat[2] <- 0
  direction <- setNames(rep(1, 3), names(omics))

  expect_warning(
    suppressMessages(ORBIT_P(omics, direction, rho = 0)),
    "clamped"
  )
  res <- suppressMessages(suppressWarnings(
    ORBIT_P(omics, direction, rho = 0)
  ))
  expect_true(all(is.finite(res$P) | is.na(res$P)))
})

test_that("NA in stat is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- sim$omics_list
  omics[[1]]$stat[1] <- NA
  direction <- setNames(rep(1, 3), names(omics))

  expect_error(
    ORBIT_P(omics, direction),
    "NA"
  )
})


# ---------------------------------------------------------------------
# F. rho boundary handling (scalar mode)
# ---------------------------------------------------------------------
test_that("Negative supplied rho is coerced to 0 with a warning", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  expect_warning(
    suppressMessages(ORBIT_P(sim$omics_list, direction, rho = -0.1)),
    "negative"
  )
  res_neg <- suppressMessages(suppressWarnings(
    ORBIT_P(sim$omics_list, direction, rho = -0.1)
  ))
  res_zero <- suppressMessages(ORBIT_P(sim$omics_list, direction, rho = 0))
  expect_equal(res_neg$P, res_zero$P)
})

test_that("rho = 1 is allowed but emits a Laplace warning", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  expect_warning(
    suppressMessages(ORBIT_P(sim$omics_list, direction, rho = 1)),
    "Laplace"
  )
  res <- suppressMessages(suppressWarnings(
    ORBIT_P(sim$omics_list, direction, rho = 1)
  ))
  expect_true(all(res$P >= 0 & res$P <= 1, na.rm = TRUE))
})

test_that("rho > 1 is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = 1.5),
    "must not exceed 1"
  )
})


# ---------------------------------------------------------------------
# G. Input validation (errors)
# ---------------------------------------------------------------------
test_that("Non-list / unnamed omics_list is rejected", {
  expect_error(
    ORBIT_P(matrix(1, 3, 3), direction = c(a = 1)),
    "named list"
  )
  unnamed <- list(
    data.frame(feature = "f1", sign = 1, stat = 0.5,
               stringsAsFactors = FALSE),
    data.frame(feature = "f1", sign = 1, stat = 0.5,
               stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_P(unnamed, direction = c(1, 1)),
    "named list"
  )
})

test_that("Non-data-frame omic element is rejected", {
  omics <- list(a = matrix(1, 3, 3),
                b = data.frame(feature = "f1", sign = 1, stat = 0.5,
                               stringsAsFactors = FALSE))
  expect_error(
    ORBIT_P(omics, c(a = 1, b = 1)),
    "not a data frame"
  )
})

test_that("Missing required column is rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2"),
                   stat = c(0.1, 0.2),  # no `sign`
                   stringsAsFactors = FALSE),
    b = data.frame(feature = c("f1", "f2"), sign = c(1, -1),
                   stat = c(0.1, 0.2),
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_P(omics, c(a = 1, b = 1)),
    "missing column"
  )
})

test_that("Bad sign values are rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 0, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   stringsAsFactors = FALSE),
    b = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 1, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_P(omics, c(a = 1, b = 1)),
    "must be -1 or \\+1"
  )
})

test_that("Duplicate feature names within an omic are rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f1", "f2"),
                   sign    = c(1, 1, 1),
                   stat    = c(0.1, 0.2, 0.3),
                   stringsAsFactors = FALSE),
    b = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 1, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_P(omics, c(a = 1, b = 1)),
    "duplicate"
  )
})

test_that("Unnamed direction is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  expect_error(
    ORBIT_P(sim$omics_list, c(1, 1, 1)),
    "named numeric"
  )
})

test_that("direction names must match omic names", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  expect_error(
    ORBIT_P(sim$omics_list, c(wrong = 1, names = 1, here = 1)),
    "must match"
  )
})

test_that("Non-+/-1 or NA direction is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  oc <- names(sim$omics_list)

  expect_error(
    ORBIT_P(sim$omics_list, setNames(c(1, 0, 1), oc)),
    "must contain only -1 or \\+1"
  )
  expect_error(
    ORBIT_P(sim$omics_list, setNames(c(1, NA, 1), oc)),
    "must contain only -1 or \\+1"
  )
})


# ---------------------------------------------------------------------
# H. Feature membership across omics
# ---------------------------------------------------------------------
test_that("Features missing from some omics get reduced N", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  omics <- sim$omics_list
  feat  <- omics[[1]]$feature

  omics[[1]] <- omics[[1]][!omics[[1]]$feature %in% feat[1:5], , drop = FALSE]
  omics[[2]] <- omics[[2]][!omics[[2]]$feature %in% feat[10:12], , drop = FALSE]

  direction <- setNames(rep(1, 3), names(omics))
  res <- suppressMessages(ORBIT_P(omics, direction, rho = 0))

  m <- match(feat, res$Feature)
  expect_equal(res$N[m[1:5]],   rep(2L, 5))
  expect_equal(res$N[m[10:12]], rep(2L, 3))
  expect_equal(res$N[m[20:30]], rep(3L, 11))
})


# ---------------------------------------------------------------------
# I. direction order-independence
# ---------------------------------------------------------------------
test_that("direction can be supplied in any order", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  oc <- names(sim$omics_list)
  direction <- setNames(c(1, 1, -1), oc)

  res1 <- suppressMessages(
    ORBIT_P(sim$omics_list, direction,             rho = 0.2))
  res2 <- suppressMessages(
    ORBIT_P(sim$omics_list, direction[c(3, 1, 2)], rho = 0.2))
  expect_equal(res1, res2)
})


# ---------------------------------------------------------------------
# J. Matrix rho mode
# ---------------------------------------------------------------------
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
    ORBIT_P(sim$omics_list, direction, rho = rho_val))
  res_matrix <- suppressMessages(
    ORBIT_P(sim$omics_list, direction, rho = rho_mat))

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
  rho_perm <- rho_mat[c(3, 1, 2), c(3, 1, 2)]

  res1 <- suppressMessages(
    ORBIT_P(sim$omics_list, direction, rho = rho_mat))
  res2 <- suppressMessages(
    ORBIT_P(sim$omics_list, direction, rho = rho_perm))
  expect_equal(res1, res2)
})

test_that("Matrix rho gives per-feature eff_rho based on presence pattern", {
  # Features 1:50 only in omics 1,2; features 51:100 in all 3.
  feat <- paste0("f", 1:100); oc <- paste0("o", 1:3)
  set.seed(1)
  stat_mat <- matrix(runif(100 * 3), 100, 3)
  direction <- setNames(rep(1, 3), oc)

  omics <- setNames(list(
    data.frame(feature = feat, sign = 1, stat = stat_mat[, 1],
               stringsAsFactors = FALSE),
    data.frame(feature = feat, sign = 1, stat = stat_mat[, 2],
               stringsAsFactors = FALSE),
    data.frame(feature = feat[51:100], sign = 1,
               stat = stat_mat[51:100, 3],
               stringsAsFactors = FALSE)
  ), oc)

  rho_mat <- matrix(c(1.0, 0.01, 0.90,
                      0.01, 1.0, 0.90,
                      0.90, 0.90, 1.0),
                    nrow = 3, dimnames = list(oc, oc))

  res <- suppressMessages(
    ORBIT_P(omics, direction, rho = rho_mat))

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
      ORBIT_P(omics, direction, rho = rho_mat)),
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
    ORBIT_P(sim$omics_list, direction, rho = bad),
    "omics"
  )

  bad2 <- matrix(0.2, 3, 4)
  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = bad2),
    "square"
  )
})

test_that("Matrix rho rejects missing / wrong names", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  oc <- names(sim$omics_list)

  bad <- matrix(0.2, 3, 3); diag(bad) <- 1
  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = bad),
    "rownames and colnames"
  )

  bad2 <- make_const_rho_mat(c("x", "y", "z"), 0.2)
  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = bad2),
    "must match"
  )
})

test_that("Matrix rho rejects non-symmetric / bad diagonal / off-diag >= 1", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  oc <- names(sim$omics_list)

  # Non-symmetric
  bad1 <- matrix(c(1.0, 0.2, 0.1,
                   0.3, 1.0, 0.4,
                   0.1, 0.4, 1.0),
                 nrow = 3, dimnames = list(oc, oc))
  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = bad1),
    "symmetric"
  )

  # Diagonal != 1
  bad2 <- make_const_rho_mat(oc, 0.2)
  bad2[1, 1] <- 0.99
  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = bad2),
    "diagonal"
  )

  # Off-diagonal = 1
  bad3 <- make_const_rho_mat(oc, 1.0)
  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = bad3),
    "strictly less than 1"
  )

  # Off-diagonal > 1
  bad4 <- make_const_rho_mat(oc, 1.5)
  expect_error(
    ORBIT_P(sim$omics_list, direction, rho = bad4),
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

  # off-diag = -0.4; inflation = 1 + 3 * (-0.4) = -0.2 <= 0
  rho_mat <- make_const_rho_mat(oc, -0.4)

  expect_warning(
    res <- suppressMessages(
      ORBIT_P(omics, direction, rho = rho_mat)),
    "non-positive"
  )
  expect_true(all(is.na(res$P)))
  expect_true(all(is.na(res$Direction)))
})
