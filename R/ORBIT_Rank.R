#' Variance-Gamma Integration of Multi-Omic Rank Statistics
#'
#' Combines per-omic signed rank statistics into a single
#' per-feature p-value via a Variance-Gamma null calibrated by the
#' average between-omic correlation. The correlation can be supplied
#' as a single scalar or as a full pairwise matrix; in the matrix
#' case the effective correlation is computed per feature from the
#' subset of omics in which that feature is observed.
#'
#' @param omics_list Named list of data frames, one per omic. Each must
#'   contain three columns whose names are given by `feature_col`,
#'   `sign_col`, `stat_col`. Different omics may report different
#'   subsets of features; the function works on the union.
#' @param direction Named numeric vector of +1/-1, names must match
#'   `names(omics_list)`. Multiplied into the score so that the user
#'   can flip an omic's sign convention without editing `sign_col`.
#' @param rho Between-omic correlation. Either:
#'   * A single number in `(-1/(k_max - 1), 1)`, where `k_max` is the
#'     largest number of omics in which any single feature is
#'     observed. The same value is used for every feature. For two
#'     omics (k_max = 2) any value in `(-1, 1)` is allowed. Negative
#'     `rho` shrinks the null variance and increases power for
#'     direction-aligned signals (e.g. testing reversal between two
#'     contrasts); positive `rho` accounts for shared technical /
#'     batch nuisance correlation.
#'   * A numeric square symmetric matrix with `rownames` and
#'     `colnames` matching `names(omics_list)` (order may differ; it
#'     is reordered internally). Diagonal must be exactly 1; off-
#'     diagonal entries must be strictly less than 1. `NA` is allowed
#'     in off-diagonal entries for pairs whose correlation could not
#'     be estimated (see "Matrix mode" below for handling).
#'   Default `0`.
#' @param stat_direction Convention for `stat_col`. `"smaller_better"`
#'   (default) treats smaller `stat` as more significant (p-value
#'   style); `"larger_better"` treats larger `stat` as more
#'   significant (`|t-stat|`, `-log10(p)`, enrichment score, etc.).
#'   Either way `stat` must be non-negative; the sign of the effect
#'   belongs in `sign_col`.
#' @param feature_col,sign_col,stat_col Column names inside each data
#'   frame. Defaults: `"feature"`, `"sign"`, `"stat"`. `feature` may
#'   be character or factor (factors are coerced); `sign` must be -1
#'   or +1 (no NA); `stat` must be numeric, non-negative, with no NA.
#' @param seed Optional seed for random tie-breaking in ranking. The
#'   global RNG state is restored on exit.
#'
#' @details
#' **Score construction.** Within each omic, features are ranked by
#' `stat` (most-significant first, controlled by `stat_direction`)
#' with random tie-breaking. The per-feature signed score is
#' `sign * direction * (-log(adj_rank))`, where
#' `adj_rank = r / (m + 1)` and `m` is the number of features in that
#' omic. Under the null `|score| ~ Exp(1)` and the signed score is
#' marginally Laplace(0, 1). A feature not observed in a given omic
#' contributes nothing to its combined statistic.
#'
#' **Combined statistic.** For each feature, `D = sum over omics of
#' score`, taken over the `k` omics in which the feature is observed.
#' `N` (returned) is this `k`. Under the equicorrelated Laplace null
#' with correlation `eff_rho`, `D` follows a symmetric Variance-Gamma
#' distribution; the correlation inflation `1 + (k - 1) * eff_rho` is
#' absorbed via moment matching (variance exactly; higher moments
#' approximate). The inflation must be strictly positive, which
#' requires `eff_rho > -1/(k - 1)`.
#'
#' **Scalar mode.** `eff_rho = rho` for every feature. The function
#' computes `k_max = max(N)` over features with `N > 0` and enforces
#' `rho > -1/(k_max - 1)` (with a small numerical tolerance);
#' violation is an error. When `rho < 0` and `k_max >= 2`, a message
#' notes that the negative value is within the valid range and
#' increases power for direction-aligned signals.
#'
#' **Matrix mode.** Each feature has its own *presence pattern* (the
#' set of omics in which it is observed). Features are grouped by
#' pattern and assigned `eff_rho` = the mean of the upper-triangular
#' entries of `rho` restricted to the omics in that pattern. Single-
#' omic features (`k = 1`) get `eff_rho = 0` as a placeholder, since
#' `(k - 1) * eff_rho = 0` regardless. Handling of `NA` entries in
#' `rho`:
#'   * If the relevant submatrix has *some* `NA` pairs but at least
#'     one finite pair, `eff_rho` is the mean of the finite pairs and
#'     a single `warning()` reports the count.
#'   * If *all* relevant pairs are `NA`, `eff_rho` falls back to `0`
#'     and the same warning reports this count separately.
#'   * If the resulting inflation `1 + (k - 1) * eff_rho` is not
#'     strictly positive (within `1e-6`), the feature is skipped:
#'     `eff_rho` stays `NA` and the feature's returned `P` is `NA`.
#'     A separate warning reports the count.
#' In matrix mode a message also reports the median and range of
#' per-feature `eff_rho` over features with `N >= 2`.
#'
#' **Tail probability.** Per-feature `P = P(|D'| >= |D|)` is computed
#' on the Variance-Gamma with shape `r_eff = k / (1 + (k-1)*eff_rho)`
#' and scale `sigma_eff = 1 + (k-1)*eff_rho`, evaluated by 32-node
#' Gauss-Laguerre quadrature on the upper tail in log space. Special
#' cases: `r_eff = 1` (Laplace) and `r_eff = 2` use closed forms.
#' Features are grouped by `(N, eff_rho)` so the quadrature setup is
#' reused.
#'
#' **Reported feature overlap.** Before scoring, the function emits
#' `message()`s with the size of the feature union, the number and
#' percentage of features shared across *all* omics, and a note
#' listing any omic pair with zero shared features (informational; it
#' does not affect this function but may be relevant to the choice of
#' `rho` if `rho` was estimated from such a pair).
#'
#' @return Data frame with columns:
#'   * `Feature`: feature identifier (union across omics).
#'   * `N`: integer, number of omics in which the feature is observed
#'     (i.e. has non-NA score). `0` is possible only if a feature
#'     appears in some omic but with no rows — in practice every
#'     returned feature has `N >= 1`.
#'   * `P`: two-sided Variance-Gamma tail probability. `NA` when
#'     `N == 0`, or (matrix mode only) when the feature's effective
#'     inflation is non-positive.
#'   * `Direction`: integer, `1L` if `D > 0`, `-1L` if `D < 0`, `NA`
#'     if `D == 0` or `P` is `NA`.
#'
#'@examples
#' # 500 features, 80% null, 20% true signal (10% up + 10% down)
#' sim <- ORBIT_simulate_signal(
#'   n = 500, K = 3,
#'   pi_null = 0.8, pi_up = 0.1, pi_down = 0.1,
#'   mu_signal = 2,
#'   rho_null = 0, rho_signal = 0.5,
#'   seed = 1
#' )
#' direction <- setNames(rep(1, 3), names(sim$omics_list))
#'
#' # Independence assumption: rho = 0 (matches rho_null here)
#' res <- ORBIT_Rank(sim$omics_list, direction, rho = 0, seed = 1)
#' head(res)
#'
#' # Calibration sanity check against the simulated truth
#' fdr   <- p.adjust(res$P, method = "BH")
#' truth <- sim$is_signal[res$Feature]
#' table(predicted = fdr < 0.05, truth = truth)
#' @export
ORBIT_Rank <- function(omics_list,
                       direction,
                       rho = 0,
                       stat_direction = c("smaller_better", "larger_better"),
                       feature_col = "feature",
                       sign_col    = "sign",
                       stat_col    = "stat",
                       seed = NULL) {

  # ---- 0. Mode + column-name parameters ----
  stat_direction <- match.arg(stat_direction)

  check_col <- function(x, nm) {
    if (!is.character(x) || length(x) != 1 || is.na(x) || x == "") {
      stop("`", nm, "` must be a single non-empty character string.")
    }
  }
  check_col(feature_col, "feature_col")
  check_col(sign_col,    "sign_col")
  check_col(stat_col,    "stat_col")
  if (anyDuplicated(c(feature_col, sign_col, stat_col))) {
    stop("`feature_col`, `sign_col`, `stat_col` must be three distinct names.")
  }

  sig_rank <- function(t) {
    if (stat_direction == "larger_better") {
      rank(-t, ties.method = "random")
    } else {
      rank( t, ties.method = "random")
    }
  }

  # ---- 1. omics_list structural checks ----
  if (!is.list(omics_list) || is.null(names(omics_list)) ||
      any(names(omics_list) == "") || anyDuplicated(names(omics_list))) {
    stop("`omics_list` must be a named list with unique, non-empty omic names.")
  }
  omic_names <- names(omics_list)
  n_omic <- length(omics_list)
  required_cols <- c(feature_col, sign_col, stat_col)

  # ---- 2. Per-omic content checks; standardize columns to feature/sign/stat ----
  for (j in seq_len(n_omic)) {
    df_j <- omics_list[[j]]
    nm <- omic_names[j]
    if (!is.data.frame(df_j)) {
      stop("Omic '", nm, "' is not a data frame.")
    }
    miss <- setdiff(required_cols, colnames(df_j))
    if (length(miss)) {
      stop("Omic '", nm, "' is missing column(s): ",
           paste(miss, collapse = ", "), ".")
    }
    f <- df_j[[feature_col]]
    s <- df_j[[sign_col]]
    t <- df_j[[stat_col]]

    if (is.factor(f)) f <- as.character(f)
    if (!is.character(f)) {
      stop("Omic '", nm, "': `", feature_col, "` must be character or factor.")
    }
    if (anyNA(f) || any(f == "")) {
      stop("Omic '", nm, "': `", feature_col, "` has missing or empty value(s).")
    }
    if (anyDuplicated(f)) {
      stop("Omic '", nm, "': duplicate values in `", feature_col, "`.")
    }
    if (!is.numeric(s)) {
      stop("Omic '", nm, "': `", sign_col, "` must be numeric.")
    }
    if (anyNA(s)) {
      stop("Omic '", nm, "': `", sign_col, "` has ", sum(is.na(s)),
           " missing value(s).")
    }
    if (any(!s %in% c(-1, 1)))
      stop("Omic '", nm, "': `", sign_col, "` must be -1 or +1.")
    if (!is.numeric(t)) {
      stop("Omic '", nm, "': `", stat_col, "` must be numeric.")
    }
    if (anyNA(t)) {
      stop("Omic '", nm, "': `", stat_col, "` has ", sum(is.na(t)),
           " missing value(s).")
    }
    neg <- t < 0
    if (any(neg)) {
      stop("Omic '", nm, "': `", stat_col, "` has ", sum(neg),
           " negative value(s). `", stat_col, "` must be non-negative; ",
           "put the sign of the effect in `", sign_col, "`.")
    }
    omics_list[[j]] <- data.frame(feature = f, sign = s, stat = t,
                                  stringsAsFactors = FALSE)
  }

  # ---- 3. rho: accept scalar OR a named symmetric matrix ----
  is_matrix_rho <- is.matrix(rho)

  if (is_matrix_rho) {
    if (!is.numeric(rho))
      stop("`rho` matrix must be numeric.")
    if (nrow(rho) != ncol(rho))
      stop("`rho` matrix must be square.")
    if (nrow(rho) != n_omic)
      stop(sprintf("`rho` matrix is %dx%d but `omics_list` has %d omics.",
                   nrow(rho), ncol(rho), n_omic))
    rn <- rownames(rho); cn <- colnames(rho)
    if (is.null(rn) || is.null(cn))
      stop("`rho` matrix must have rownames and colnames matching `omics_list` names.")
    if (!setequal(rn, omic_names) || !setequal(cn, omic_names))
      stop("`rho` matrix rownames/colnames must match `omics_list` names.")
    rho <- rho[omic_names, omic_names, drop = FALSE]

    ut_vals <- rho[upper.tri(rho)]
    lt_vals <- t(rho)[upper.tri(rho)]
    na_mismatch <- is.na(ut_vals) != is.na(lt_vals)
    val_mismatch <- !is.na(ut_vals) & !is.na(lt_vals) &
      abs(ut_vals - lt_vals) > 1e-8
    if (any(na_mismatch) || any(val_mismatch))
      stop("`rho` matrix must be symmetric.")

    diag_rho <- diag(rho)
    if (any(is.na(diag_rho)) || any(abs(diag_rho - 1) > 1e-8))
      stop("`rho` matrix diagonal must be exactly 1.")

    off_finite <- ut_vals[is.finite(ut_vals)]
    if (length(off_finite) && any(off_finite >= 1 - 1e-12))
      stop("Off-diagonal `rho` values must be strictly less than 1.")

  } else {
    if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho))
      stop("`rho` must be a single finite number or a numeric matrix.")
    if (rho > 1) stop("`rho` must not exceed 1, got ", rho, ".")
    if (rho == 1) {
      stop(sprintf("`rho = %.4f` is not allowed; must be strictly less than 1. ", rho),
           "At rho = 1 the omics are perfectly correlated under the null, ",
           "so combining them provides no information beyond a single omic.")
    }
  }

  # ---- 4. direction ----
  if (!is.numeric(direction) || is.null(names(direction))) {
    stop("`direction` must be a named numeric vector.")
  }
  if (!setequal(names(direction), omic_names)) {
    stop("Names of `direction` must match names of `omics_list`. ",
         "Missing: ", paste(setdiff(omic_names, names(direction)),
                            collapse = ", "), ". ",
         "Extra: ",   paste(setdiff(names(direction), omic_names),
                            collapse = ", "), ".")
  }
  direction <- direction[omic_names]
  if (anyNA(direction) || !all(direction %in% c(-1, 1))) {
    stop("`direction` must contain only -1 or +1, no NA.")
  }

  # ---- 5. Seed handling ----
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
      stop("`seed` must be a single finite number or NULL.")
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (is.null(old_seed)) {
        suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  # ---- 6. Feature-overlap report ----
  feature_lists <- lapply(omics_list, `[[`, "feature")
  all_features  <- unique(unlist(feature_lists, use.names = FALSE))
  n_total <- length(all_features)
  shared <- if (n_omic >= 1) Reduce(intersect, feature_lists) else character(0)
  n_shared <- length(shared)
  pct_shared <- if (n_total > 0) 100 * n_shared / n_total else 0

  message(sprintf("Union of features across %d omics: %d",
                  n_omic, n_total))
  message(sprintf("Shared across ALL omics: %d (%.1f%% of union)",
                  n_shared, pct_shared))

  if (n_omic >= 2) {
    empty_pairs <- character(0)
    for (a in seq_len(n_omic - 1)) {
      for (b in seq.int(a + 1, n_omic)) {
        if (length(intersect(feature_lists[[a]], feature_lists[[b]])) == 0) {
          empty_pairs <- c(empty_pairs,
                           paste0(omic_names[a], " <-> ", omic_names[b]))
        }
      }
    }
    if (length(empty_pairs)) {
      message("Note: Omic pair(s) with no shared features (does not affect ORBIT_Rank): ",
              paste(empty_pairs, collapse = "; "), ". ",
              "If overlap was expected, please verify that all omics use a common ",
              "feature ID space (e.g., consistent ID types like gene symbol vs Ensembl, ",
              "capitalization, organism).")
    }
  }

  # ---- 7. Signed-score matrix ----
  signed_score_mat <- matrix(NA_real_, n_total, n_omic,
                             dimnames = list(all_features, omic_names))
  for (j in seq_len(n_omic)) {
    df_j <- omics_list[[j]]
    m_j <- nrow(df_j)
    if (m_j == 0) next
    r_j <- sig_rank(df_j$stat)
    adj_rank_j <- r_j / (m_j + 1)
    score_j <- df_j$sign * direction[j] * (-log(adj_rank_j))
    signed_score_mat[df_j$feature, j] <- score_j
  }
  N_vec <- rowSums(!is.na(signed_score_mat))
  D_vec <- rowSums(signed_score_mat, na.rm = TRUE)

  # ---- 7b. Effective rho per feature; validity check ----
  eff_rho_vec <- rep(NA_real_, n_total)

  if (is_matrix_rho) {
    presence_mat <- !is.na(signed_score_mat)
    pattern_keys <- apply(presence_mat, 1L,
                          function(r) paste(which(r), collapse = ","))
    n_partial_na <- 0L
    n_all_na     <- 0L
    n_bad_inf    <- 0L
    eps          <- 1e-6
    for (pat in unique(pattern_keys)) {
      if (!nzchar(pat)) next
      idx <- which(pattern_keys == pat)
      omics_set <- as.integer(strsplit(pat, ",", fixed = TRUE)[[1]])
      k <- length(omics_set)
      if (k == 1L) {
        eff_rho_vec[idx] <- 0
        next
      }
      sub_rho   <- rho[omics_set, omics_set, drop = FALSE]
      pair_vals <- sub_rho[upper.tri(sub_rho)]
      n_pairs   <- length(pair_vals)
      n_na      <- sum(is.na(pair_vals))
      if (n_na == 0L) {
        mr <- mean(pair_vals)
      } else if (n_na < n_pairs) {
        mr <- mean(pair_vals, na.rm = TRUE)
        n_partial_na <- n_partial_na + length(idx)
      } else {
        mr <- 0                                  # all pairs NA -> default
        n_all_na <- n_all_na + length(idx)
      }
      if (1 + (k - 1) * mr <= eps) {
        n_bad_inf <- n_bad_inf + length(idx)
        next                                     # eff_rho stays NA -> P = NA
      }
      eff_rho_vec[idx] <- mr
    }
    if (n_partial_na > 0L || n_all_na > 0L) {
      warning(sprintf(
        "%d feature(s) had some NA pairs in `rho` matrix (effective rho computed from available pairs); %d feature(s) had ALL pairs NA (defaulted to rho = 0).",
        n_partial_na, n_all_na))
    }
    if (n_bad_inf > 0L) {
      warning(sprintf(
        "%d feature(s) skipped due to non-positive inflation (effective rho too negative for their k). Their P is NA.",
        n_bad_inf))
    }
    finite_eff <- eff_rho_vec[is.finite(eff_rho_vec) & N_vec >= 2L]
    if (length(finite_eff)) {
      message(sprintf(
        "Matrix-mode rho: per-feature effective rho over %d features with N>=2: median=%.4f, range=[%.4f, %.4f]",
        length(finite_eff),
        stats::median(finite_eff), min(finite_eff), max(finite_eff)))
    }
  } else {
    k_pos <- N_vec[N_vec > 0]
    k_max <- if (length(k_pos)) max(k_pos) else 1L
    if (k_max >= 2) {
      rho_lb <- -1 / (k_max - 1)
      eps <- 1e-6
      if (rho <= rho_lb + eps) {
        stop(sprintf(
          "`rho` = %.4f is at or below the valid lower bound %.4f for k_max = %d. The equicorrelated Laplace null requires rho > -1/(k_max - 1), i.e. inflation > 0.",
          rho, rho_lb, k_max))
      }
      if (n_omic >= 3) {
        message(sprintf(
          "Using scalar rho = %.4f for every feature, regardless of which omic subset it is observed in. With %d omics this is a uniform-correlation assumption; if pairwise correlations across omic pairs differ substantially, the same scalar applied to all features can bias the null calibration. Pass `rho` as a named symmetric matrix (e.g. `ORBIT_cor()$rho_mat`) for per-feature effective rho.",
          rho, n_omic))
      }
      if (rho < 0) {
        message(sprintf(
          "Using negative rho = %.4f (valid range for these data: (%.4f, 1)). Null variance is reduced; this increases power for direction-aligned signals.",
          rho, rho_lb))
      }
    }
    eff_rho_vec[N_vec > 0] <- rho
  }

  # ---- 8. Variance-Gamma tail probability ----
  GL <- statmod::gauss.quad(32, kind = "laguerre")
  GL_NODES <- GL$nodes
  GL_LOG_W <- log(GL$weights)

  log_besselK_safe <- function(z, nu) {
    if (!is.finite(z) || z <= 0) return(NA_real_)
    val <- suppressWarnings(besselK(z, nu, expon.scaled = TRUE))
    if (is.finite(val) && val > 0) return(log(val) - z)
    NA_real_
  }
  log_vg_density_local <- function(ax, r, sigma) {
    if (!is.finite(ax) || ax <= 0) return(NA_real_)
    nu <- r - 0.5
    z  <- ax / sigma
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
      u <- GL_NODES[i]
      x <- abs_d + sigma * u
      lg <- log_vg_density_local(x, r, sigma)
      if (!is.finite(lg)) return(NA_real_)
      lg + u + GL_LOG_W[i]
    }, numeric(1))
    if (any(is.na(log_vals))) return(NA_real_)
    m <- max(log_vals)
    if (!is.finite(m)) return(NA_real_)
    log(2) + log(sigma) + m + log(sum(exp(log_vals - m)))
  }

  # ---- 9. P per feature, grouped by (N, eff_rho) ----
  logP_vec <- rep(NA_real_, n_total)
  valid_idx <- which(N_vec > 0L & !is.na(eff_rho_vec))
  if (length(valid_idx)) {
    grp_key <- paste(N_vec[valid_idx],
                     format(eff_rho_vec[valid_idx], digits = 15),
                     sep = "_")
    for (key in unique(grp_key)) {
      sel       <- valid_idx[grp_key == key]
      k_val     <- N_vec[sel[1]]
      eff_val   <- eff_rho_vec[sel[1]]
      inflation <- 1 + (k_val - 1) * eff_val
      r_eff     <- k_val / inflation
      sigma_eff <- inflation
      logP_vec[sel] <- vapply(
        abs(D_vec[sel]),
        function(d) get_log_tail_prob(d, r_eff, sigma_eff),
        numeric(1)
      )
    }
  }
  P_vec <- exp(logP_vec)

  # ---- 10. direction: integer -1/+1, NA when D = 0 or P unavailable ----
  dir_vec <- as.integer(sign(D_vec))
  dir_vec[dir_vec == 0] <- NA_integer_
  dir_vec[is.na(P_vec)] <- NA_integer_

  # ---- 11. Output ----
  data.frame(
    Feature   = all_features,
    N         = as.integer(N_vec),
    P         = P_vec,
    Direction = dir_vec,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
