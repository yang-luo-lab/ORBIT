#' Variance-Gamma Integration of Multi-Omic P-Values
#'
#' Combines per-feature P-values across omics using `-log(P)`
#' directly as the signed evidence score (no within-omic ranking).
#' Use this when you trust P-value magnitudes and consider them
#' comparable across omics. Otherwise feed the same P-values to
#' `ORBIT_Rank()`, which ranks within each omic and is invariant
#' to magnitude differences. Supports per-feature effective
#' between-omic correlation when `rho` is supplied as a matrix or
#' auto-estimated.
#'
#' @param omics_list Named list of data frames, one per omic. Each
#'   must contain columns named by `feature_col`, `sign_col`,
#'   `stat_col`. Different omics may report different subsets of
#'   features; the function works on the union. `stat_col` is
#'   interpreted as a P-value in `[0, 1]`.
#' @param direction Named numeric vector of +1/-1, names must match
#'   `names(omics_list)`.
#' @param rho Between-omic correlation. One of:
#'   * `NULL` (default): auto-estimate a pairwise correlation matrix
#'     using the same signed rank-based scoring as `ORBIT_cor()`,
#'     applied to all features. Skipping signal filtering is
#'     acceptable here because P-value combination via `-log(P)` is
#'     dominated by large signal magnitudes, making the combined
#'     statistic relatively robust to slight `rho` misspecification —
#'     more so than rank-based combination, where each feature's
#'     contribution is capped at `log(m+1)`. This mirrors Empirical
#'     Brown's Method, which also estimates correlation without
#'     explicit filtering. For a cleaner estimate that excludes
#'     significant features, run `ORBIT_cor()` and pass its
#'     `$rho_mat` here.
#'   * A single number in `[0, 1]`. The same value is applied to
#'     every feature. Negative user-supplied values are coerced to
#'     `0` with a warning; values above `1` are an error; `rho = 1`
#'     triggers a warning (VG null degenerates to Laplace).
#'   * A numeric square symmetric matrix with `rownames` / `colnames`
#'     matching `names(omics_list)`. Validated and used exactly as
#'     in `ORBIT_Rank()`'s matrix mode: diagonal must be `1`,
#'     off-diagonal must be strictly less than `1`, `NA` off-diagonal
#'     entries allowed.
#'   When `rho` is a matrix (user-supplied or auto-estimated), each
#'   feature gets its own effective `rho` equal to the mean of the
#'   upper-triangular entries of `rho` restricted to the omics in
#'   which that feature is observed (its "presence pattern").
#' @param feature_col,sign_col,stat_col Column names. Defaults
#'   `"feature"`, `"sign"`, `"stat"`. `feature` may be character or
#'   factor; `sign` must be -1 or +1 (no NA); `stat` must be a
#'   P-value in `[0, 1]` (no NA). P-values equal to exactly `0` are
#'   clamped to `.Machine$double.xmin` with a single aggregated
#'   warning.
#' @param seed Optional seed for random tie-breaking in the rank
#'   transform used by the auto-estimator. The global RNG state is
#'   restored on exit.
#'
#' @details
#' **Combined statistic.** Per-feature signed score in omic `j` is
#' `sign * direction[j] * (-log(P))`. Under the null,
#' `|score| ~ Exp(1)`. The combined statistic `D = sum of scores` over
#' the `k = N` observed omics follows, under the equicorrelated
#' Laplace null with effective correlation `eff_rho`, a symmetric
#' Variance-Gamma with shape `r_eff = k / (1 + (k - 1) * eff_rho)`
#' and scale `sigma_eff = 1 + (k - 1) * eff_rho`. The tail probability
#' uses 32-node Gauss-Laguerre quadrature, with closed forms at
#' `r_eff = 1` (Laplace) and `r_eff = 2`. Features are grouped by
#' `(N, eff_rho)` so the quadrature setup is reused.
#'
#' **Per-feature effective rho (matrix mode).** Each feature uses the
#' mean of upper-triangular `rho` entries restricted to its observed
#' omics. Single-omic features get `eff_rho = 0`. NA entries in `rho`
#' use available pairs (with `rho = 0` as fallback when all are NA).
#' Features whose `1 + (k - 1) * eff_rho` is non-positive are skipped
#' (`P = NA`).

#' @return Data frame with columns:
#'   * `Feature`: feature identifier (union across omics).
#'   * `N`: integer, number of omics in which the feature is observed.
#'   * `P`: two-sided Variance-Gamma tail probability. `NA` when
#'     `N == 0` or when the feature is skipped due to non-positive
#'     effective inflation.
#'   * `Direction`: integer, `1L` if `D > 0`, `-1L` if `D < 0`, `NA`
#'     if `D == 0` or `P` is `NA`.
#'   The `rho` actually used (scalar or matrix) is attached as
#'   `attr(result, "rho")` for inspection.
#'
#' @importFrom stats cor median quantile
#' @examples
#' # 500 features, 80% null, 20% true signal, with non-zero null
#' # between-omic correlation so auto-estimate has something to recover.
#' sim <- ORBIT_simulate_signal(
#'   n = 500, K = 3,
#'   pi_null = 0.8, pi_up = 0.1, pi_down = 0.1,
#'   mu_signal = 2,
#'   rho_null = 0.3, rho_signal = 0.5,
#'   seed = 1
#' )
#' direction <- setNames(rep(1, 3), names(sim$omics_list))
#'
#' # rho not supplied: auto-estimated as a matrix and attached
#' res <- ORBIT_P(sim$omics_list, direction)
#' attr(res, "rho")
#' head(res)
#'
#' # Calibration sanity check against the simulated truth
#' fdr   <- p.adjust(res$P, method = "BH")
#' truth <- sim$is_signal[res$Feature]
#' table(predicted = fdr < 0.05, truth = truth)
#' @export
ORBIT_P <- function(omics_list,
                    direction,
                    rho = NULL,
                    feature_col = "feature",
                    sign_col    = "sign",
                    stat_col    = "stat",
                    seed        = NULL) {

  # ---- 0. Column-name parameter checks ----
  check_col <- function(x, nm) {
    if (!is.character(x) || length(x) != 1 || is.na(x) || x == "")
      stop("`", nm, "` must be a single non-empty character string.")
  }
  check_col(feature_col, "feature_col")
  check_col(sign_col,    "sign_col")
  check_col(stat_col,    "stat_col")
  if (anyDuplicated(c(feature_col, sign_col, stat_col)))
    stop("`feature_col`, `sign_col`, `stat_col` must be three distinct names.")

  # ---- 1. omics_list structural checks ----
  if (!is.list(omics_list) || is.null(names(omics_list)) ||
      any(names(omics_list) == "") || anyDuplicated(names(omics_list)))
    stop("`omics_list` must be a named list with unique, non-empty omic names.")
  omic_names <- names(omics_list); n_omic <- length(omics_list)
  required_cols <- c(feature_col, sign_col, stat_col)

  # ---- 2. Per-omic content checks; standardize columns ----
  total_zero_clamped <- 0L
  for (j in seq_len(n_omic)) {
    df_j <- omics_list[[j]]; nm <- omic_names[j]
    if (!is.data.frame(df_j)) stop("Omic '", nm, "' is not a data frame.")
    miss <- setdiff(required_cols, colnames(df_j))
    if (length(miss))
      stop("Omic '", nm, "' is missing column(s): ",
           paste(miss, collapse = ", "), ".")
    f <- df_j[[feature_col]]; s <- df_j[[sign_col]]; t <- df_j[[stat_col]]
    if (is.factor(f)) f <- as.character(f)
    if (!is.character(f))
      stop("Omic '", nm, "': `", feature_col, "` must be character or factor.")
    if (anyNA(f) || any(f == ""))
      stop("Omic '", nm, "': `", feature_col, "` has missing or empty value(s).")
    if (anyDuplicated(f))
      stop("Omic '", nm, "': duplicate values in `", feature_col, "`.")
    if (!is.numeric(s)) stop("Omic '", nm, "': `", sign_col, "` must be numeric.")
    if (anyNA(s))      stop("Omic '", nm, "': `", sign_col, "` has NA.")
    if (any(!s %in% c(-1, 1)))
      stop("Omic '", nm, "': `", sign_col, "` must be -1 or +1.")
    if (!is.numeric(t)) stop("Omic '", nm, "': `", stat_col, "` must be numeric.")
    if (anyNA(t))       stop("Omic '", nm, "': `", stat_col, "` has NA.")
    if (any(t < 0))
      stop("Omic '", nm, "': `", stat_col,
           "` (P-values) must be in [0, 1]; negative values found.")
    if (any(t > 1))
      stop("Omic '", nm, "': `", stat_col,
           "` (P-values) must not exceed 1. Use `ORBIT_Rank()` for rank stats.")
    n_zero <- sum(t == 0)
    if (n_zero > 0) {
      t[t == 0] <- .Machine$double.xmin
      total_zero_clamped <- total_zero_clamped + n_zero
    }
    omics_list[[j]] <- data.frame(feature = f, sign = s, stat = t,
                                  stringsAsFactors = FALSE)
  }
  if (total_zero_clamped > 0)
    warning(total_zero_clamped, " P-value(s) equal to 0 were clamped to ",
            format(.Machine$double.xmin, digits = 3), ".")

  # ---- 3. rho input dispatch ----
  auto_estimate <- is.null(rho)
  is_matrix_rho <- !auto_estimate && is.matrix(rho)

  if (is_matrix_rho) {
    if (!is.numeric(rho))         stop("`rho` matrix must be numeric.")
    if (nrow(rho) != ncol(rho))   stop("`rho` matrix must be square.")
    if (nrow(rho) != n_omic)
      stop(sprintf("`rho` matrix is %dx%d but `omics_list` has %d omics.",
                   nrow(rho), ncol(rho), n_omic))
    rn <- rownames(rho); cn <- colnames(rho)
    if (is.null(rn) || is.null(cn))
      stop("`rho` matrix must have rownames and colnames matching omic names.")
    if (!setequal(rn, omic_names) || !setequal(cn, omic_names))
      stop("`rho` matrix rownames/colnames must match `omics_list` names.")
    rho <- rho[omic_names, omic_names, drop = FALSE]
    ut_vals <- rho[upper.tri(rho)]
    lt_vals <- t(rho)[upper.tri(rho)]
    if (any(is.na(ut_vals) != is.na(lt_vals)) ||
        any(!is.na(ut_vals) & !is.na(lt_vals) &
            abs(ut_vals - lt_vals) > 1e-8))
      stop("`rho` matrix must be symmetric.")
    diag_rho <- diag(rho)
    if (any(is.na(diag_rho)) || any(abs(diag_rho - 1) > 1e-8))
      stop("`rho` matrix diagonal must be exactly 1.")
    off_finite <- ut_vals[is.finite(ut_vals)]
    if (length(off_finite) && any(off_finite >= 1 - 1e-12))
      stop("Off-diagonal `rho` values must be strictly less than 1 (max = ",
           max(off_finite), ").")
  } else if (!auto_estimate) {
    if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho))
      stop("`rho` must be a single finite number, a numeric matrix, or NULL.")
    if (rho > 1) stop("`rho` must not exceed 1, got ", rho, ".")
    if (rho < 0) {
      warning("`rho` was negative (", rho, "); coerced to 0.")
      rho <- 0
    }
    if (rho == 1)
      warning("`rho = 1`: Variance-Gamma null degenerates to Laplace.")
  }

  # ---- 4. direction ----
  if (!is.numeric(direction) || is.null(names(direction)))
    stop("`direction` must be a named numeric vector.")
  if (!setequal(names(direction), omic_names))
    stop("Names of `direction` must match names of `omics_list`.")
  direction <- direction[omic_names]
  if (anyNA(direction) || !all(direction %in% c(-1, 1)))
    stop("`direction` must contain only -1 or +1, no NA.")

  # ---- 5. Seed handling ----
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed))
      stop("`seed` must be a single finite number or NULL.")
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old_seed))
        suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
      else assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  # ---- 6. Feature-overlap report ----
  feature_lists <- lapply(omics_list, `[[`, "feature")
  all_features  <- unique(unlist(feature_lists, use.names = FALSE))
  n_total <- length(all_features)
  shared  <- if (n_omic >= 1) Reduce(intersect, feature_lists) else character(0)
  pct_shared <- if (n_total > 0) 100 * length(shared) / n_total else 0
  message(sprintf("Union of features across %d omics: %d", n_omic, n_total))
  message(sprintf("Shared across ALL omics: %d (%.1f%% of union)",
                  length(shared), pct_shared))
  if (n_omic >= 2) {
    empty_pairs <- character(0)
    for (a in seq_len(n_omic - 1)) for (b in seq.int(a + 1, n_omic)) {
      if (length(intersect(feature_lists[[a]], feature_lists[[b]])) == 0)
        empty_pairs <- c(empty_pairs,
                         paste0(omic_names[a], " <-> ", omic_names[b]))
    }
    if (length(empty_pairs))
      warning("Omic pair(s) with NO shared features: ",
              paste(empty_pairs, collapse = "; "),
              ". If overlap was expected, please verify that all omics use a ",
              "common feature ID space (e.g., consistent ID types, ",
              "capitalization, organism).")
  }

  # ---- 7. Signed-score matrix using -log(P) directly (for D) + wide stat/sign ----
  signed_score_mat <- matrix(NA_real_, n_total, n_omic,
                             dimnames = list(all_features, omic_names))
  stat_wide <- matrix(NA_real_, n_total, n_omic,
                      dimnames = list(all_features, omic_names))
  sign_wide <- matrix(NA_real_, n_total, n_omic,
                      dimnames = list(all_features, omic_names))
  for (j in seq_len(n_omic)) {
    df_j <- omics_list[[j]]; if (!nrow(df_j)) next
    idx <- match(df_j$feature, all_features)
    stat_wide[idx, j] <- df_j$stat
    sign_wide[idx, j] <- df_j$sign
    signed_score_mat[idx, j] <- df_j$sign * direction[j] * (-log(df_j$stat))
  }
  has_feat <- !is.na(signed_score_mat)
  N_vec <- rowSums(has_feat)
  D_vec <- rowSums(signed_score_mat, na.rm = TRUE)

  # ---- 8. Auto-estimate rho matrix: per-pair rerank within overlap, Pearson ----
  if (auto_estimate) {
    if (n_omic < 2) {
      rho <- 0; is_matrix_rho <- FALSE
      message("Single-omic input; setting rho = 0.")
    } else {
      rho_est <- diag(1, n_omic)
      dimnames(rho_est) <- list(omic_names, omic_names)

      for (ii in seq_len(n_omic - 1)) for (jj in seq.int(ii + 1, n_omic)) {
        overlap_ij <- has_feat[, ii] & has_feat[, jj]
        n_ov <- sum(overlap_ij)
        if (n_ov < 2L) {
          rho_est[ii, jj] <- rho_est[jj, ii] <- NA_real_
          next
        }
        stat_i <- stat_wide[overlap_ij, ii]; sign_i <- sign_wide[overlap_ij, ii]
        stat_j <- stat_wide[overlap_ij, jj]; sign_j <- sign_wide[overlap_ij, jj]
        r_i <- rank(stat_i, ties.method = "random")
        r_j <- rank(stat_j, ties.method = "random")
        adj_i <- r_i / (n_ov + 1)
        adj_j <- r_j / (n_ov + 1)
        score_i <- sign_i * direction[ii] * (-log(adj_i))
        score_j <- sign_j * direction[jj] * (-log(adj_j))

        r_t <- suppressWarnings(stats::cor(score_i, score_j, method = "pearson"))
        rho_est[ii, jj] <- rho_est[jj, ii] <-
          if (is.finite(r_t)) r_t else NA_real_
      }

      ut_vals      <- rho_est[upper.tri(rho_est)]
      finite_pairs <- ut_vals[is.finite(ut_vals)]
      n_na_pairs   <- sum(is.na(ut_vals))

      if (length(finite_pairs) == 0L) {
        rho <- 0; is_matrix_rho <- FALSE
        warning("Could not estimate any pairwise correlation; falling back to rho = 0.")
      } else {
        cap_thr <- 1 - 1e-12
        ut_mask <- upper.tri(rho_est)
        bad_high <- which(ut_vals > cap_thr & !is.na(ut_vals))
        if (length(bad_high)) {
          rho_est[ut_mask][bad_high] <- cap_thr
          rho_est <- (rho_est + t(rho_est)) / 2
          diag(rho_est) <- 1
          warning(sprintf("%d auto-estimated off-diagonal correlation(s) >= 1 capped.",
                          length(bad_high)))
        }
        rho <- rho_est; is_matrix_rho <- TRUE
        message(sprintf(
          "Auto-estimated rho matrix (per-pair signed rank-based, no signal filtering): %d/%d finite pairs, median=%.4f, range=[%.4f, %.4f].",
          length(finite_pairs), length(ut_vals),
          stats::median(finite_pairs), min(finite_pairs), max(finite_pairs)))
        if (n_na_pairs > 0L)
          message(sprintf("  %d pair(s) had NA correlation.", n_na_pairs))
      }
    }
  }

  # ---- 9. Per-feature effective rho ----
  eff_rho_vec <- rep(NA_real_, n_total)

  if (is_matrix_rho) {
    pattern_keys <- apply(has_feat, 1L,
                          function(r) paste(which(r), collapse = ","))
    n_partial_na <- 0L; n_all_na <- 0L; n_bad_inf <- 0L
    eps <- 1e-6
    for (pat in unique(pattern_keys)) {
      if (!nzchar(pat)) next
      idx <- which(pattern_keys == pat)
      omics_set <- as.integer(strsplit(pat, ",", fixed = TRUE)[[1]])
      k <- length(omics_set)
      if (k == 1L) { eff_rho_vec[idx] <- 0; next }
      sub_rho   <- rho[omics_set, omics_set, drop = FALSE]
      pair_vals <- sub_rho[upper.tri(sub_rho)]
      n_pairs   <- length(pair_vals); n_na <- sum(is.na(pair_vals))
      if (n_na == 0L) {
        mr <- mean(pair_vals)
      } else if (n_na < n_pairs) {
        mr <- mean(pair_vals, na.rm = TRUE)
        n_partial_na <- n_partial_na + length(idx)
      } else {
        mr <- 0
        n_all_na <- n_all_na + length(idx)
      }
      if (1 + (k - 1) * mr <= eps) {
        n_bad_inf <- n_bad_inf + length(idx); next
      }
      eff_rho_vec[idx] <- mr
    }
    if (n_partial_na > 0L || n_all_na > 0L)
      warning(sprintf(
        "%d feature(s) had some NA pairs (used available); %d feature(s) had ALL pairs NA (defaulted to rho = 0).",
        n_partial_na, n_all_na))
    if (n_bad_inf > 0L)
      warning(sprintf(
        "%d feature(s) skipped due to non-positive effective inflation. Their P is NA.",
        n_bad_inf))
    finite_eff <- eff_rho_vec[is.finite(eff_rho_vec) & N_vec >= 2L]
    if (length(finite_eff))
      message(sprintf(
        "Per-feature effective rho over %d features with N>=2: median=%.4f, range=[%.4f, %.4f]",
        length(finite_eff), stats::median(finite_eff),
        min(finite_eff), max(finite_eff)))
  } else {
    eff_rho_vec[N_vec > 0] <- rho
  }

  # ---- 10. Variance-Gamma tail probability ----
  GL <- statmod::gauss.quad(32, kind = "laguerre")
  GL_NODES <- GL$nodes; GL_LOG_W <- log(GL$weights)
  log_besselK_safe <- function(z, nu) {
    if (!is.finite(z) || z <= 0) return(NA_real_)
    val <- suppressWarnings(besselK(z, nu, expon.scaled = TRUE))
    if (is.finite(val) && val > 0) return(log(val) - z)
    NA_real_
  }
  log_vg_density_local <- function(ax, r, sigma) {
    if (!is.finite(ax) || ax <= 0) return(NA_real_)
    nu <- r - 0.5; z <- ax / sigma
    log_pref  <- -log(sigma) - 0.5 * log(pi) - lgamma(r)
    log_power <- nu * (log(ax) - log(2 * sigma))
    logK <- log_besselK_safe(z, nu)
    if (!is.finite(logK)) return(NA_real_)
    log_pref + log_power + logK
  }
  get_log_tail_prob <- function(abs_d, r, sigma) {
    if (!is.finite(abs_d) || r <= 0 || sigma <= 0) return(NA_real_)
    if (abs_d <= 1e-10) return(0)
    if (abs(r - 1) < 1e-3) return(-abs_d / sigma)
    if (abs(r - 2) < 1e-3) {
      z <- abs_d / sigma
      return(log(0.5) + log(z + 2) - z)
    }
    log_vals <- vapply(seq_along(GL_NODES), function(i) {
      u <- GL_NODES[i]; x <- abs_d + sigma * u
      lg <- log_vg_density_local(x, r, sigma)
      if (!is.finite(lg)) return(NA_real_)
      lg + u + GL_LOG_W[i]
    }, numeric(1))
    if (any(is.na(log_vals))) return(NA_real_)
    m <- max(log_vals)
    if (!is.finite(m)) return(NA_real_)
    log(2) + log(sigma) + m + log(sum(exp(log_vals - m)))
  }

  # ---- 11. P per feature, grouped by (N, eff_rho) ----
  logP_vec <- rep(NA_real_, n_total)
  valid_idx <- which(N_vec > 0L & !is.na(eff_rho_vec))
  if (length(valid_idx)) {
    grp_key <- paste(N_vec[valid_idx],
                     format(eff_rho_vec[valid_idx], digits = 15),
                     sep = "_")
    for (key in unique(grp_key)) {
      sel       <- valid_idx[grp_key == key]
      k_val     <- N_vec[sel[1]]; eff_val <- eff_rho_vec[sel[1]]
      inflation <- 1 + (k_val - 1) * eff_val
      r_eff     <- k_val / inflation; sigma_eff <- inflation
      logP_vec[sel] <- vapply(
        abs(D_vec[sel]),
        function(d) get_log_tail_prob(d, r_eff, sigma_eff),
        numeric(1))
    }
  }
  P_vec <- exp(logP_vec)

  # ---- 12. Direction: integer -1/+1, NA when D = 0 or P unavailable ----
  dir_vec <- as.integer(sign(D_vec))
  dir_vec[dir_vec == 0]  <- NA_integer_
  dir_vec[is.na(P_vec)]  <- NA_integer_

  # ---- 13. Output ----
  result <- data.frame(
    Feature   = all_features,
    N         = as.integer(N_vec),
    P         = P_vec,
    Direction = dir_vec,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(result, "rho") <- rho   # scalar or matrix actually used
  result
}
