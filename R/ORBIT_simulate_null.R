#' Simulate Multi-Omic Null Data with Specified Between-Omic Correlation
#'
#' Generates a synthetic multi-omic dataset under the null hypothesis:
#' all features are non-significant, but P-values across omics are
#' allowed to share a common pairwise correlation `rho`. The output
#' is in the `omics_list` format expected by [ORBIT_P()],
#' [ORBIT_Rank()], and [ORBIT_cor()], so it can be fed directly to
#' any of those without further conversion. Optional MCAR
#' missingness can be injected by dropping features from selected
#' omics' data frames.
#'
#' Internally, a length-`K` standard multivariate normal vector is
#' drawn for each feature using a Cholesky decomposition of the
#' exchangeable correlation matrix (1 on the diagonal, `rho` off the
#' diagonal). Each component is converted to a two-sided normal
#' P-value, with the corresponding sign retained. If `missing_rate`
#' is positive, a fraction of features in each target omic is
#' removed from that omic's data frame entirely (completely at random
#' within each targeted omic). No external packages are required for
#' the simulation.
#'
#' @param n Number of features to simulate. Must be a positive
#'   integer.
#' @param K Number of omics. Must be at least 2.
#' @param rho Average pairwise correlation between omics, in the
#'   open interval `(-1 / (K - 1), 1)`. Defaults to `0`
#'   (independent omics).
#' @param missing_rate Fraction of features dropped from each target
#'   omic. Must be in `[0, 1)`. Defaults to `0` (no missingness).
#'   Missingness is independent across target omics (MCAR).
#' @param missing_omics Which omics to apply missingness to. Either
#'   the string `"all"` (every omic, default), the string `"last"`
#'   (only the last omic), or an integer vector of omic indices in
#'   `1:K`. Ignored when `missing_rate = 0`.
#' @param seed Optional integer seed for reproducibility. If `NULL`
#'   (default), the current random state is used.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`omics_list`}{Named list of length `K`, one data frame
#'       per omic, with columns `feature` (character), `sign`
#'       (`+1` / `-1`), and `stat` (two-sided P-value). Names are
#'       `omic1` ... `omicK`. Ready to pass directly to [ORBIT_P()],
#'       [ORBIT_Rank()], or [ORBIT_cor()]. Features dropped by
#'       `missing_rate` simply do not appear in the corresponding
#'       data frame.}
#'     \item{`Z_mat`}{The underlying `n` by `K` standard normal
#'       draws, with row names `feature1` ... `featureN` and column
#'       names `omic1` ... `omicK`. Entries corresponding to dropped
#'       features are `NA` so callers can recover which features
#'       were removed from each omic.}
#'     \item{`Sigma`}{The `K` by `K` exchangeable correlation matrix
#'       used for sampling.}
#'   }
#'
#' @export
#'
#' @examples
#' # Independent null, ready to feed into ORBIT_P
#' sim <- ORBIT_simulate_null(n = 200, K = 3, rho = 0, seed = 1)
#' direction <- setNames(rep(1, 3), names(sim$omics_list))
#' res <- ORBIT_P(sim$omics_list, direction)
#' head(res)
#'
#' # Correlated null
#' sim2 <- ORBIT_simulate_null(n = 200, K = 4, rho = 0.4, seed = 1)
#' length(sim2$omics_list)
#'
#' # Half of features missing in the last omic only
#' sim3 <- ORBIT_simulate_null(n = 200, K = 2, rho = 0.3,
#'                             missing_rate = 0.5,
#'                             missing_omics = "last",
#'                             seed = 1)
#' sapply(sim3$omics_list, nrow)
ORBIT_simulate_null <- function(n = 1000,
                                K = 2,
                                rho = 0,
                                missing_rate = 0,
                                missing_omics = "all",
                                seed = NULL) {

  # ----------------------------------------------------------
  # 1. Argument checking
  # ----------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != as.integer(n)) {
    stop("`n` must be a positive integer.")
  }
  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != as.integer(K)) {
    stop("`K` must be an integer >= 2.")
  }
  if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho)) {
    stop("`rho` must be a single finite number.")
  }
  rho_lower <- -1 / (K - 1)
  if (rho <= rho_lower || rho >= 1) {
    stop(sprintf(
      "`rho` must lie in the open interval (%.4f, 1) for K = %d.",
      rho_lower, K
    ))
  }
  if (!is.numeric(missing_rate) || length(missing_rate) != 1 ||
      !is.finite(missing_rate) || missing_rate < 0 ||
      missing_rate >= 1) {
    stop("`missing_rate` must be a single number in [0, 1).")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("`seed` must be a single number or NULL.")
    }
    set.seed(seed)
  }

  n <- as.integer(n)
  K <- as.integer(K)

  # Resolve `missing_omics` into a vector of column indices
  if (is.character(missing_omics) && length(missing_omics) == 1) {
    missing_omics <- match.arg(missing_omics, c("all", "last"))
    target_cols <- if (missing_omics == "all") seq_len(K) else K
  } else if (is.numeric(missing_omics)) {
    if (any(!is.finite(missing_omics)) ||
        any(missing_omics != as.integer(missing_omics)) ||
        any(missing_omics < 1) || any(missing_omics > K)) {
      stop("Numeric `missing_omics` must contain integer indices ",
           "in 1:K.")
    }
    target_cols <- as.integer(unique(missing_omics))
  } else {
    stop("`missing_omics` must be 'all', 'last', or an integer vector.")
  }

  # ----------------------------------------------------------
  # 2. Build exchangeable correlation matrix and its Cholesky factor
  # ----------------------------------------------------------
  Sigma <- matrix(rho, nrow = K, ncol = K)
  diag(Sigma) <- 1
  L <- chol(Sigma)   # upper-triangular: t(L) %*% L = Sigma

  # ----------------------------------------------------------
  # 3. Draw correlated standard normals
  # ----------------------------------------------------------
  Z_iid <- matrix(stats::rnorm(n * K), nrow = n, ncol = K)
  Z <- Z_iid %*% L

  # ----------------------------------------------------------
  # 4. Two-sided P-values and signs
  # ----------------------------------------------------------
  P <- 2 * stats::pnorm(-abs(Z))
  P <- pmax(P, .Machine$double.xmin)
  sign_mat <- ifelse(Z >= 0, 1, -1)

  # ----------------------------------------------------------
  # 5. Names
  # ----------------------------------------------------------
  feat_names <- paste0("feature", seq_len(n))
  omic_names <- paste0("omic",    seq_len(K))
  dn <- list(feat_names, omic_names)
  dimnames(P)        <- dn
  dimnames(sign_mat) <- dn
  dimnames(Z)        <- dn
  dimnames(Sigma)    <- list(omic_names, omic_names)

  # ----------------------------------------------------------
  # 6. MCAR missingness: choose features to drop per target omic
  # ----------------------------------------------------------
  drop_mask <- matrix(FALSE, nrow = n, ncol = K, dimnames = dn)
  if (missing_rate > 0 && length(target_cols) > 0) {
    n_miss <- round(n * missing_rate)
    if (n_miss > 0) {
      for (j in target_cols) {
        miss_idx <- sample.int(n, size = n_miss, replace = FALSE)
        drop_mask[miss_idx, j] <- TRUE
      }
    }
  }
  # Reflect drops in Z_mat (for the diagnostic field only)
  Z[drop_mask] <- NA_real_

  # ----------------------------------------------------------
  # 7. Build per-omic data frames; dropped features simply absent
  # ----------------------------------------------------------
  omics_list <- setNames(lapply(seq_len(K), function(j) {
    keep <- !drop_mask[, j]
    data.frame(
      feature = feat_names[keep],
      sign    = sign_mat[keep, j],
      stat    = P[keep, j],
      stringsAsFactors = FALSE
    )
  }), omic_names)

  list(
    omics_list = omics_list,
    Z_mat      = Z,
    Sigma      = Sigma
  )
}
