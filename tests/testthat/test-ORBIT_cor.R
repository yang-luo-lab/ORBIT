# ====================================================================
# Tests for ORBIT_cor
# ====================================================================

# Helper: add significance = 0 column to a null simulator's output
add_zero_sig <- function(omics_list) {
  lapply(omics_list, function(df) {
    df$significance <- 0L
    df
  })
}

# ---------------------------------------------------------------------
# A. Basic structure of the return value
# ---------------------------------------------------------------------
test_that("ORBIT_cor returns a list with the expected structure", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  res <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  )

  expect_type(res, "list")
  expect_named(res,
               c("rho", "rho_mat", "diag",
                 "n_total", "tukey_quartiles", "settings"))
  expect_length(res$rho, 1)

  expect_true(is.matrix(res$rho_mat))
  expect_equal(dim(res$rho_mat), c(3, 3))
  expect_equal(rownames(res$rho_mat), names(omics))
  expect_equal(colnames(res$rho_mat), names(omics))

  expect_s3_class(res$diag, "data.frame")
  expect_named(res$diag,
               c("omic_i", "omic_j", "n_overlap", "n_bg",
                 "bg_frac", "rho", "tukey_dev", "rho_loo"))
  expect_equal(nrow(res$diag), 3L)   # 3 choose 2 = 3 pairs

  expect_type(res$n_total, "integer")
  expect_gt(res$n_total, 0L)

  expect_named(res$tukey_quartiles, c("Q1", "Q3", "IQR"))

  expect_type(res$settings, "list")
  expect_named(res$settings, c("min_bg", "bg_frac_tol"))
})

test_that("rho is finite and within [-1, 1]", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  res <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  )
  expect_true(is.finite(res$rho))
  expect_gte(res$rho, -1); expect_lte(res$rho, 1)
})

test_that("rho_mat diagonal is 1", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  res <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  )
  expect_true(all(abs(diag(res$rho_mat) - 1) < 1e-8))
})

test_that("rho_mat is symmetric", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  res <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  )
  expect_true(isSymmetric(res$rho_mat))
})

test_that("rho equals mean of upper-triangle finite cors", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  res <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  )
  ut_vals <- res$rho_mat[upper.tri(res$rho_mat)]
  expected_rho <- mean(ut_vals[is.finite(ut_vals)])
  expect_equal(res$rho, expected_rho, tolerance = 1e-12)
})


# ---------------------------------------------------------------------
# B. Statistical correctness of rho estimation
# ---------------------------------------------------------------------
test_that("Independent null: estimated rho is near 0", {
  sim <- ORBIT_simulate_null(n = 3000, K = 3, rho = 0, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))
  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  expect_lt(abs(res$rho), 0.005)
})

test_that("Correlated null: estimated rho recovers the truth", {
  true_rho <- 0.4
  sim <- ORBIT_simulate_null(n = 3000, K = 3, rho = true_rho, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))
  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  expect_true(is.finite(res$rho))
  expect_lt(abs(res$rho - true_rho), 0.005)
})

test_that("Mixture data: ORBIT_cor recovers rho_null on the background", {
  sim <- ORBIT_simulate_signal(
    n = 3000, K = 3,
    pi_null = 0.95, pi_up = 0.025, pi_down = 0.025,
    mu_signal = 5,
    rho_null = 0.2, rho_signal = 0.5,
    seed = 1
  )
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(sim$omics_list, direction,
              sig_mode = "thr", sig_thr = 0.05, seed = 1)
  ))
  expect_lt(abs(res$rho - 0.2), 0.01)
})


# ---------------------------------------------------------------------
# C. Significance specification
# ---------------------------------------------------------------------
test_that("sig_mode = 'column' uses the supplied significance column", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  sig_features <- paste0("feature", 1:5)
  omics <- lapply(sim$omics_list, function(df) {
    df$significance <- as.integer(df$feature %in% sig_features)
    df
  })

  res <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  )
  N <- nrow(sim$omics_list[[1]])
  expect_true(all(res$diag$n_bg == N - 5L))
})

test_that("stat_direction = 'larger_better' selects top by magnitude", {
  set.seed(7)
  feat <- paste0("feature", 1:60)
  ocol <- paste0("omic", 1:3)
  rank_stat <- abs(matrix(rnorm(60 * 3), 60, 3))
  sign_mat  <- matrix(sample(c(-1, 1), 60 * 3, replace = TRUE), 60, 3)
  direction <- setNames(rep(1, 3), ocol)

  omics <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat,
               sign    = sign_mat[, j],
               stat    = rank_stat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction,
              sig_mode = "thr", sig_thr = 0.10,
              stat_direction = "larger_better", seed = 1)
  ))
  expect_true(is.finite(res$rho))
  expect_true(all(res$diag$n_bg >= 60 - 12))
})

test_that("stat_direction flips which end is treated as significant", {
  set.seed(42)
  feat <- paste0("f", 1:200); ocol <- paste0("o", 1:2)
  stat_mat <- matrix(runif(400), 200, 2)
  direction <- setNames(c(1, 1), ocol)

  omics <- setNames(lapply(seq_along(ocol), function(j) {
    data.frame(feature = feat, sign = 1, stat = stat_mat[, j],
               stringsAsFactors = FALSE)
  }), ocol)

  r_small <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "thr", sig_thr = 0.1,
              stat_direction = "smaller_better", seed = 1)))
  r_large <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "thr", sig_thr = 0.1,
              stat_direction = "larger_better", seed = 1)))
  expect_true(is.finite(r_small$rho))
  expect_true(is.finite(r_large$rho))
  expect_true(all(r_small$diag$n_bg > 0))
  expect_true(all(r_large$diag$n_bg > 0))
})


# ---------------------------------------------------------------------
# D. Input validation (error paths)
# ---------------------------------------------------------------------
test_that("Non-list / unnamed omics_list is rejected", {
  expect_error(
    ORBIT_cor(matrix(1, 3, 3),
              direction = c(a = 1),
              sig_mode = "thr", sig_thr = 0.05),
    "named list"
  )
  unnamed <- list(
    data.frame(feature = "f1", sign = 1, stat = 0.5,
               stringsAsFactors = FALSE),
    data.frame(feature = "f1", sign = 1, stat = 0.5,
               stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_cor(unnamed, direction = c(1, 1),
              sig_mode = "thr", sig_thr = 0.05),
    "named list"
  )
})

test_that("Fewer than 2 omics is rejected", {
  omics <- list(a = data.frame(feature = "f1", sign = 1, stat = 0.5,
                               stringsAsFactors = FALSE))
  expect_error(
    ORBIT_cor(omics, c(a = 1), sig_mode = "thr", sig_thr = 0.5),
    "at least 2 omics"
  )
})

test_that("Non-data-frame omic element is rejected", {
  omics <- list(a = matrix(1, 3, 3),
                b = data.frame(feature = "f1", sign = 1, stat = 0.5,
                               stringsAsFactors = FALSE))
  expect_error(
    ORBIT_cor(omics, c(a = 1, b = 1), sig_mode = "thr", sig_thr = 0.5),
    "not a data frame"
  )
})

test_that("Missing required column is rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2"),
                   stat = c(0.1, 0.2),
                   stringsAsFactors = FALSE),
    b = data.frame(feature = c("f1", "f2"), sign = c(1, -1),
                   stat = c(0.1, 0.2), significance = 0L,
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_cor(omics, c(a = 1, b = 1),
              sig_mode = "column"),
    "missing column"
  )
})

test_that("Negative stat values are rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, -1, 1),
                   stat    = c(0.1, -0.2, 0.3),
                   significance = 0L,
                   stringsAsFactors = FALSE),
    b = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 1, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   significance = 0L,
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_cor(omics, c(a = 1, b = 1), sig_mode = "column"),
    "negative values"
  )
})

test_that("Bad sign values are rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 0, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   significance = 0L,
                   stringsAsFactors = FALSE),
    b = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 1, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   significance = 0L,
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_cor(omics, c(a = 1, b = 1), sig_mode = "column"),
    "must be -1 or \\+1"
  )
})

test_that("Duplicate feature names within an omic are rejected", {
  omics <- list(
    a = data.frame(feature = c("f1", "f1", "f2"),
                   sign    = c(1, 1, 1),
                   stat    = c(0.1, 0.2, 0.3),
                   significance = 0L,
                   stringsAsFactors = FALSE),
    b = data.frame(feature = c("f1", "f2", "f3"),
                   sign    = c(1, 1, -1),
                   stat    = c(0.1, 0.2, 0.3),
                   significance = 0L,
                   stringsAsFactors = FALSE)
  )
  expect_error(
    ORBIT_cor(omics, c(a = 1, b = 1), sig_mode = "column"),
    "duplicate"
  )
})

test_that("Unnamed direction is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  expect_error(
    ORBIT_cor(omics, c(1, 1, 1), sig_mode = "column"),
    "named numeric"
  )
})

test_that("direction names must match omic names", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  expect_error(
    ORBIT_cor(omics, c(wrong = 1, names = 1, here = 1),
              sig_mode = "column"),
    "must match"
  )
})

test_that("Invalid sig_thr is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  expect_error(
    ORBIT_cor(sim$omics_list, direction, sig_mode = "thr", sig_thr = 1.5),
    "in \\(0, 1\\)"
  )
  expect_error(
    ORBIT_cor(sim$omics_list, direction, sig_mode = "thr", sig_thr = 0),
    "in \\(0, 1\\)"
  )
})

test_that("sig_mode = 'thr' without sig_thr is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  expect_error(
    ORBIT_cor(sim$omics_list, direction, sig_mode = "thr"),
    "sig_thr"
  )
})

test_that("Invalid min_bg is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))
  expect_error(
    ORBIT_cor(omics, direction, sig_mode = "column", min_bg = 1),
    "min_bg.*>= 2"
  )
})

test_that("Invalid bg_frac_tol is rejected", {
  sim <- ORBIT_simulate_null(n = 50, K = 3, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))
  expect_error(
    ORBIT_cor(omics, direction, sig_mode = "column", bg_frac_tol = -0.1),
    "bg_frac_tol.*\\[0, 1\\]"
  )
  expect_error(
    ORBIT_cor(omics, direction, sig_mode = "column", bg_frac_tol = 1.5),
    "bg_frac_tol.*\\[0, 1\\]"
  )
})


# ---------------------------------------------------------------------
# E. direction order-independence
# ---------------------------------------------------------------------
test_that("direction can be supplied in any order", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  oc  <- names(omics)
  direction <- setNames(c(1, 1, -1), oc)

  res1 <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1))
  res2 <- suppressMessages(
    ORBIT_cor(omics, direction[c(3, 1, 2)], sig_mode = "column", seed = 1))
  expect_equal(res1$rho, res2$rho)
  expect_equal(res1$rho_mat, res2$rho_mat)
})


# ---------------------------------------------------------------------
# F. Edge cases & reproducibility
# ---------------------------------------------------------------------
test_that("Too few background features yields NA rho with warning", {
  sim <- ORBIT_simulate_null(n = 20, K = 3, seed = 1)
  omics <- lapply(sim$omics_list, function(df) {
    df$significance <- 1L
    df
  })
  direction <- setNames(rep(1, 3), names(omics))

  expect_warning(
    res <- suppressMessages(
      ORBIT_cor(omics, direction, sig_mode = "column")),
    "bg <"
  )
  expect_true(is.na(res$rho))
  # All 3 pairs should have n_bg = 0 and rho_mat off-diag = NA
  expect_true(all(res$diag$n_bg == 0L))
  ut_vals <- res$rho_mat[upper.tri(res$rho_mat)]
  expect_true(all(is.na(ut_vals)))
})

test_that("Same seed yields reproducible results", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  r1 <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 42))
  r2 <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 42))
  expect_identical(r1$rho, r2$rho)
  expect_identical(r1$rho_mat, r2$rho_mat)
  expect_identical(r1$diag, r2$diag)
})

test_that("Global .Random.seed is restored after a seeded call", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  set.seed(123)
  before <- .Random.seed
  suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 999))
  after <- .Random.seed
  expect_identical(before, after)
})


# ---------------------------------------------------------------------
# G. min_bg parameter
# ---------------------------------------------------------------------
test_that("min_bg controls per-pair background filtering", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))
  r1 <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column",
              min_bg = 2L, seed = 1))
  expect_true(is.finite(r1$rho))
  expect_true(all(is.finite(r1$rho_mat[upper.tri(r1$rho_mat)])))
  r2 <- suppressWarnings(suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column",
              min_bg = 1e6, seed = 1)))
  expect_true(is.na(r2$rho))
  expect_true(all(is.na(r2$rho_mat[upper.tri(r2$rho_mat)])))
})

test_that("settings echoes the supplied parameters", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 3), names(omics))

  res <- suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column",
              min_bg = 50L, bg_frac_tol = 0.3, seed = 1))
  expect_equal(res$settings$min_bg, 50L)
  expect_equal(res$settings$bg_frac_tol, 0.3)
})


# ---------------------------------------------------------------------
# H. bg_frac_tol parameter
# ---------------------------------------------------------------------
test_that("bg_frac_tol = 0 silences the low-fraction warning", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  omics <- lapply(sim$omics_list, function(df) {
    df$significance <- 0L
    df$significance[1:70] <- 1L
    df
  })
  direction <- setNames(rep(1, 3), names(omics))
  expect_no_warning(suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column",
              bg_frac_tol = 0, seed = 1)
  ))
})

test_that("bg_frac_tol = 1 fires warning when any bg is below overlap", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, seed = 1)
  omics <- lapply(sim$omics_list, function(df) {
    df$significance <- 0L
    df$significance[1:10] <- 1L
    df
  })
  direction <- setNames(rep(1, 3), names(omics))

  expect_warning(suppressMessages(
    ORBIT_cor(omics, direction, sig_mode = "column",
              bg_frac_tol = 1, seed = 1)
  ), "bg fraction")
})


# ---------------------------------------------------------------------
# I. diag table structure & sorting
# ---------------------------------------------------------------------
test_that("diag rows are sorted by |tukey_dev| descending", {
  sim <- ORBIT_simulate_null(n = 500, K = 5, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 5), names(omics))

  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  td <- abs(res$diag$tukey_dev)
  td[is.na(td)] <- 0
  expect_equal(td, sort(td, decreasing = TRUE))
})

test_that("diag row count equals number of unordered pairs", {
  for (K in 2:5) {
    sim <- ORBIT_simulate_null(n = 200, K = K, seed = 1)
    omics <- add_zero_sig(sim$omics_list)
    direction <- setNames(rep(1, K), names(omics))
    res <- suppressMessages(suppressWarnings(
      ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
    ))
    expect_equal(nrow(res$diag), K * (K - 1) / 2,
                 info = paste("K =", K))
  }
})

test_that("diag bg_frac equals n_bg / n_overlap", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, seed = 1)
  omics <- lapply(sim$omics_list, function(df) {
    df$significance <- 0L
    df$significance[1:20] <- 1L
    df
  })
  direction <- setNames(rep(1, 3), names(omics))

  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  expect_equal(res$diag$bg_frac,
               res$diag$n_bg / pmax(res$diag$n_overlap, 1),
               tolerance = 1e-12)
})

test_that("diag rho values are within [-1, 1] or NA", {
  sim <- ORBIT_simulate_null(n = 200, K = 4, rho = 0.3, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 4), names(omics))

  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  rho_vec <- res$diag$rho
  finite_rho <- rho_vec[is.finite(rho_vec)]
  expect_true(all(finite_rho >= -1 & finite_rho <= 1))
})


# ---------------------------------------------------------------------
# J. tukey_dev and rho_loo computation
# ---------------------------------------------------------------------
test_that("tukey_dev and rho_loo are NA when fewer than 4 valid pairs", {
  # Only 2 omics -> 1 pair total, n_valid = 1
  sim <- ORBIT_simulate_null(n = 200, K = 2, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 2), names(omics))

  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  expect_true(all(is.na(res$diag$tukey_dev)))
  expect_true(all(is.na(res$diag$rho_loo)))    # need >= 2 valid pairs
  expect_true(all(is.na(res$tukey_quartiles))) # NA when < 4 pairs
})

test_that("rho_loo for a row equals mean of OTHER valid pair cors", {
  sim <- ORBIT_simulate_null(n = 500, K = 4, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 4), names(omics))

  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  finite_idx <- which(is.finite(res$diag$rho))
  for (i in finite_idx) {
    other <- res$diag$rho[finite_idx[finite_idx != i]]
    expect_equal(res$diag$rho_loo[i], mean(other),
                 tolerance = 1e-12,
                 info = paste("row", i))
  }
})

test_that("tukey_dev is 0 for pairs inside the IQR box", {
  sim <- ORBIT_simulate_null(n = 1000, K = 5, rho = 0, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 5), names(omics))

  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  expect_true(any(res$diag$tukey_dev == 0, na.rm = TRUE))
})

test_that("tukey_quartiles match the quartiles of finite diag rho", {
  sim <- ORBIT_simulate_null(n = 500, K = 5, rho = 0.2, seed = 1)
  omics <- add_zero_sig(sim$omics_list)
  direction <- setNames(rep(1, 5), names(omics))
  res <- suppressMessages(suppressWarnings(
    ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
  ))
  finite_rho <- res$diag$rho[is.finite(res$diag$rho)]
  if (length(finite_rho) >= 4) {
    qq <- stats::quantile(finite_rho, c(0.25, 0.75), names = FALSE)
    expect_equal(res$tukey_quartiles[["Q1"]],  qq[1], tolerance = 1e-12)
    expect_equal(res$tukey_quartiles[["Q3"]],  qq[2], tolerance = 1e-12)
    expect_equal(res$tukey_quartiles[["IQR"]], qq[2] - qq[1],
                 tolerance = 1e-12)
  }
})
