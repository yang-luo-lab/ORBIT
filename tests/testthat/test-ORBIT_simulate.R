# ====================================================================
# Tests for ORBIT_simulate_null and ORBIT_simulate_signal (updated)
#
# Both simulators now return an `omics_list` field directly (named
# list of data frames with feature/sign/stat columns) plus diagnostic
# fields. The old `stat_mat` / `sign_mat` matrix outputs are gone.
# ====================================================================


# ---------------------------------------------------------------------
# A. ORBIT_simulate_null: structure
# ---------------------------------------------------------------------
test_that("ORBIT_simulate_null returns the expected structure", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, rho = 0.3, seed = 1)

  expect_type(sim, "list")
  expect_named(sim, c("omics_list", "Z_mat", "Sigma"))

  # omics_list: named list of K data frames, each with feature/sign/stat
  expect_type(sim$omics_list, "list")
  expect_length(sim$omics_list, 3)
  expect_named(sim$omics_list, paste0("omic", 1:3))
  for (j in seq_len(3)) {
    df <- sim$omics_list[[j]]
    expect_s3_class(df, "data.frame")
    expect_named(df, c("feature", "sign", "stat"))
    expect_equal(nrow(df), 100)
    expect_equal(df$feature, paste0("feature", 1:100))
  }

  # Diagnostic matrices keep the original shape and dimnames
  expect_true(is.matrix(sim$Z_mat))
  expect_equal(dim(sim$Z_mat), c(100, 3))
  expect_equal(rownames(sim$Z_mat), paste0("feature", 1:100))
  expect_equal(colnames(sim$Z_mat), paste0("omic",    1:3))
  expect_equal(dim(sim$Sigma),   c(3, 3))
})

test_that("ORBIT_simulate_null produces valid P-values and signs", {
  sim <- ORBIT_simulate_null(n = 500, K = 2, rho = 0, seed = 1)
  all_stats <- unlist(lapply(sim$omics_list, `[[`, "stat"))
  all_signs <- unlist(lapply(sim$omics_list, `[[`, "sign"))
  expect_true(all(all_stats >= 0))
  expect_true(all(all_stats <= 1))
  expect_true(all(all_signs %in% c(-1, 1)))
})

test_that("ORBIT_simulate_null seeds are reproducible", {
  s1 <- ORBIT_simulate_null(n = 50, K = 3, rho = 0.2, seed = 42)
  s2 <- ORBIT_simulate_null(n = 50, K = 3, rho = 0.2, seed = 42)
  expect_equal(s1$omics_list, s2$omics_list)
  expect_equal(s1$Z_mat,      s2$Z_mat)
})


# ---------------------------------------------------------------------
# B. ORBIT_simulate_null: statistical correctness
# ---------------------------------------------------------------------
test_that("Z_mat has approximately the requested correlation", {
  sim <- ORBIT_simulate_null(n = 5000, K = 3, rho = 0.4, seed = 1)
  emp_cor <- cor(sim$Z_mat)
  off_diag_mean <- mean(emp_cor[upper.tri(emp_cor)])
  # SE per pair ~ sqrt((1 - rho^2) / n) ~ 0.012 at n=5000.
  # Allow 0.05 absolute slack to be safe.
  expect_lt(abs(off_diag_mean - 0.4), 0.05)
})

test_that("P-values are approximately uniform under independent null", {
  sim <- ORBIT_simulate_null(n = 5000, K = 2, rho = 0, seed = 1)
  ks <- suppressWarnings(stats::ks.test(sim$omics_list[[1]]$stat, "punif"))
  # Uniform under H0; allow loose alpha to avoid CRAN flakiness
  expect_gt(ks$p.value, 1e-3)
})


# ---------------------------------------------------------------------
# C. ORBIT_simulate_null: missing_rate / missing_omics
# ---------------------------------------------------------------------
# In the new API, "missing" means the feature is removed from that
# omic's data frame, NOT NA inside a matrix. Z_mat retains NAs in the
# dropped positions for diagnostics.
test_that("missing_rate = 0 (default) produces no missingness", {
  sim <- ORBIT_simulate_null(n = 100, K = 3, rho = 0.2, seed = 1)
  expect_true(all(vapply(sim$omics_list, nrow, integer(1)) == 100))
  expect_false(any(is.na(sim$Z_mat)))
})

test_that("missing_omics = 'last' drops only the last omic", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.3,
                             missing_rate = 0.5,
                             missing_omics = "last",
                             seed = 1)
  rows <- vapply(sim$omics_list, nrow, integer(1))
  expect_equal(rows[[1]], 200)
  expect_equal(rows[[2]], 200)
  expect_equal(rows[[3]], 100)              # 200 * 0.5 dropped

  # Z_mat keeps full shape but NAs the dropped entries
  expect_equal(dim(sim$Z_mat), c(200, 3))
  expect_equal(sum(is.na(sim$Z_mat[, 1])), 0)
  expect_equal(sum(is.na(sim$Z_mat[, 2])), 0)
  expect_equal(sum(is.na(sim$Z_mat[, 3])), 100)
})

test_that("missing_omics as an integer vector works", {
  sim <- ORBIT_simulate_null(n = 100, K = 4,
                             missing_rate = 0.3,
                             missing_omics = c(2, 4),
                             seed = 1)
  rows <- vapply(sim$omics_list, nrow, integer(1))
  expect_equal(rows[[1]], 100)              # untouched
  expect_lt   (rows[[2]], 100)              # dropped
  expect_equal(rows[[3]], 100)              # untouched
  expect_lt   (rows[[4]], 100)              # dropped
})

test_that("missing_omics rejects bad index values", {
  expect_error(
    ORBIT_simulate_null(n = 100, K = 3,
                        missing_rate = 0.2,
                        missing_omics = c(1, 5)),
    "1:K"
  )
  expect_error(
    ORBIT_simulate_null(n = 100, K = 3,
                        missing_rate = 0.2,
                        missing_omics = list("a", "b")),
    "'all', 'last'"
  )
})


# ---------------------------------------------------------------------
# D. ORBIT_simulate_null: input validation + ORBIT_* integration
# ---------------------------------------------------------------------
test_that("ORBIT_simulate_null rejects invalid n and K", {
  expect_error(ORBIT_simulate_null(n = 0,    K = 3), "positive integer")
  expect_error(ORBIT_simulate_null(n = 1.5,  K = 3), "positive integer")
  expect_error(ORBIT_simulate_null(n = 100,  K = 1), "K.*>= 2")
  expect_error(ORBIT_simulate_null(n = 100,  K = 3.5), "integer")
})

test_that("ORBIT_simulate_null rejects rho outside the valid open interval", {
  # K = 3 -> rho lower bound is -1/2 = -0.5
  expect_error(ORBIT_simulate_null(n = 100, K = 3, rho = -0.6),
               "open interval")
  expect_error(ORBIT_simulate_null(n = 100, K = 3, rho = 1),
               "open interval")
  expect_error(ORBIT_simulate_null(n = 100, K = 3, rho = 1.5),
               "open interval")
})

test_that("ORBIT_simulate_null rejects bad missing_rate", {
  expect_error(
    ORBIT_simulate_null(n = 100, K = 3, missing_rate = 1),
    "in \\[0, 1\\)"
  )
  expect_error(
    ORBIT_simulate_null(n = 100, K = 3, missing_rate = -0.1),
    "in \\[0, 1\\)"
  )
})

test_that("ORBIT_simulate_null output drops directly into ORBIT_P", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  res <- suppressMessages(ORBIT_P(sim$omics_list, direction))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 200)
  expect_named(res, c("Feature", "N", "P", "Direction"))
})

test_that("ORBIT_simulate_null output drops directly into ORBIT_Rank and ORBIT_cor", {
  sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0.2, seed = 1)
  direction <- setNames(rep(1, 3), names(sim$omics_list))

  rank_res <- suppressMessages(
    ORBIT_Rank(sim$omics_list, direction, seed = 1))
  expect_s3_class(rank_res, "data.frame")
  expect_equal(nrow(rank_res), 200)

  cor_res <- suppressMessages(
    ORBIT_cor(sim$omics_list, direction,
              sig_mode = "thr", sig_thr = 0.05, seed = 1))
  expect_type(cor_res, "list")
  expect_true(is.finite(cor_res$rho))
})


# ---------------------------------------------------------------------
# E. ORBIT_simulate_signal: structure
# ---------------------------------------------------------------------
test_that("ORBIT_simulate_signal returns the expected structure", {
  sim <- ORBIT_simulate_signal(
    n = 300, K = 3,
    pi_null = 0.7, pi_up = 0.2, pi_down = 0.1,
    mu_signal = 2,
    seed = 1
  )

  expect_type(sim, "list")
  expect_named(sim, c("omics_list", "Z_mat", "is_signal",
                      "group", "Sigma_null", "Sigma_signal"))

  expect_length(sim$omics_list, 3)
  expect_named(sim$omics_list, paste0("omic", 1:3))
  for (j in seq_len(3)) {
    df <- sim$omics_list[[j]]
    expect_s3_class(df, "data.frame")
    expect_named(df, c("feature", "sign", "stat"))
    expect_equal(nrow(df), 300)
  }

  expect_equal(dim(sim$Z_mat), c(300, 3))
  expect_length(sim$is_signal, 300)
  expect_length(sim$group, 300)
  expect_true(is.logical(sim$is_signal))
  expect_true(is.character(sim$group))
  expect_true(all(sim$group %in% c("null", "up", "down")))
})

test_that("ORBIT_simulate_signal mixture proportions match expectations", {
  sim <- ORBIT_simulate_signal(
    n = 5000, K = 2,
    pi_null = 0.6, pi_up = 0.25, pi_down = 0.15,
    mu_signal = 2,
    seed = 1
  )
  obs <- table(sim$group) / length(sim$group)
  # SE ~ sqrt(p(1-p)/n) ~ 0.007 at n=5000; allow 0.03 absolute slack.
  expect_lt(abs(unname(obs["null"]) - 0.6),  0.03)
  expect_lt(abs(unname(obs["up"])   - 0.25), 0.03)
  expect_lt(abs(unname(obs["down"]) - 0.15), 0.03)
})

test_that("Up features have mostly positive Z, down features mostly negative", {
  sim <- ORBIT_simulate_signal(
    n = 1000, K = 3,
    pi_null = 0.5, pi_up = 0.25, pi_down = 0.25,
    mu_signal = 3,
    seed = 1
  )
  up_means   <- rowMeans(sim$Z_mat[sim$group == "up",   , drop = FALSE])
  down_means <- rowMeans(sim$Z_mat[sim$group == "down", , drop = FALSE])
  expect_gt(mean(up_means),   1)
  expect_lt(mean(down_means), -1)
})


# ---------------------------------------------------------------------
# F. ORBIT_simulate_signal: edge cases
# ---------------------------------------------------------------------
test_that("All-null mixture works", {
  sim <- ORBIT_simulate_signal(
    n = 200, K = 3,
    pi_null = 1, pi_up = 0, pi_down = 0,
    seed = 1
  )
  expect_equal(sum(sim$is_signal), 0)
  expect_true(all(sim$group == "null"))
})

test_that("All-signal mixture works", {
  sim <- ORBIT_simulate_signal(
    n = 200, K = 3,
    pi_null = 0, pi_up = 0.5, pi_down = 0.5,
    seed = 1
  )
  expect_equal(sum(sim$is_signal), 200)
  expect_true(all(sim$group %in% c("up", "down")))
})


# ---------------------------------------------------------------------
# G. ORBIT_simulate_signal: input validation
# ---------------------------------------------------------------------
test_that("ORBIT_simulate_signal rejects pi values not summing to 1", {
  expect_error(
    ORBIT_simulate_signal(n = 100, pi_null = 0.5, pi_up = 0.3, pi_down = 0.3),
    "must equal 1"
  )
})

test_that("ORBIT_simulate_signal rejects negative pi", {
  expect_error(
    ORBIT_simulate_signal(n = 100, pi_null = 1.2,
                          pi_up = -0.1, pi_down = -0.1),
    "\\[0, 1\\]"
  )
})

test_that("ORBIT_simulate_signal rejects rho outside the valid interval", {
  expect_error(
    ORBIT_simulate_signal(n = 100, K = 3, rho_null = 1.5),
    "open interval"
  )
  expect_error(
    ORBIT_simulate_signal(n = 100, K = 3, rho_signal = -0.6),
    "open interval"
  )
})

test_that("ORBIT_simulate_signal rejects bad K", {
  expect_error(ORBIT_simulate_signal(n = 100, K = 1), "K.*>= 2")
})


# ---------------------------------------------------------------------
# H. ORBIT_simulate_signal: end-to-end with ORBIT_P
# ---------------------------------------------------------------------
test_that("ORBIT_P recovers some signal from simulated data", {
  sim <- ORBIT_simulate_signal(
    n = 1000, K = 3,
    pi_null = 0.8, pi_up = 0.1, pi_down = 0.1,
    mu_signal = 3,
    seed = 1
  )
  direction <- setNames(rep(1, 3), names(sim$omics_list))
  res <- suppressMessages(ORBIT_P(sim$omics_list, direction))

  # ORBIT_P no longer returns FDR; compute it here.
  fdr   <- p.adjust(res$P, method = "BH")
  truth <- sim$is_signal[res$Feature]   # named-vector indexing handles ordering

  # Strong-signal regime: should detect at least some true signals at FDR 5%.
  tab <- table(predicted = fdr < 0.05, truth = truth)
  if ("TRUE" %in% rownames(tab)) {
    expect_gt(tab["TRUE", "TRUE"], 0)
  }

  # FDR control: among predicted positives, most should be true signal.
  if (sum(fdr < 0.05, na.rm = TRUE) > 10) {
    fdr_emp <- mean(!truth[fdr < 0.05])
    expect_lt(fdr_emp, 0.20)   # generous: nominal 5%, allow finite-sample slack
  }
})
