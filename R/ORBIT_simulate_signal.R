#' Simulate Multi-Omic Mixture Data with Up- and Down-Regulated Signals
#'
#' Generates a synthetic multi-omic dataset as a three-component
#' mixture of null, up-regulated, and down-regulated features.
#' Within each component, a length-`K` multivariate normal vector is
#' drawn for each feature (with a configurable between-omic
#' correlation), then converted to two-sided P-values and signs. The
#' output is in the `omics_list` format expected by [ORBIT_P()],
#' [ORBIT_Rank()], and [ORBIT_cor()], so it can be fed directly to
#' any of those without further conversion.
#'
#' For each feature, a label is drawn from `c("null", "up", "down")`
#' with probabilities `c(pi_null, pi_up, pi_down)`. Null features
#' have mean-zero normal `Z`; up-regulated features have mean
#' `+mu_signal`; down-regulated features have mean `-mu_signal`. The
#' between-omic correlation is set to `rho_null` for null features
#' and `rho_signal` for signal features (up and down combined). All
#' random draws use base R only, so no external packages are
#' required for the simulation.
#'
#' @param n Number of features to simulate. Must be a positive
#'   integer.
#' @param K Number of omics. Must be at least 2.
#' @param pi_null Mixture proportion of null features. Must be in
#'   `[0, 1]` and `pi_null + pi_up + pi_down` must equal 1.
#' @param pi_up Mixture proportion of up-regulated features.
#' @param pi_down Mixture proportion of down-regulated features.
#' @param mu_signal Mean of the multivariate normal `Z` for up-
#'   regulated features (down-regulated features use `-mu_signal`).
#'   Larger values produce stronger signals.
#' @param rho_null Average pairwise correlation between omics among
#'   null features. Must lie in `(-1 / (K - 1), 1)`. Defaults to `0`.
#' @param rho_signal Average pairwise correlation between omics
#'   among signal features (up and down combined). Must lie in
#'   `(-1 / (K - 1), 1)`. Defaults to `0`.
#' @param seed Optional integer seed for reproducibility. If `NULL`
#'   (default), the current random state is used.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`omics_list`}{Named list of length `K`, one data frame
#'       per omic, with columns `feature` (character), `sign`
#'       (`+1` / `-1`), and `stat` (two-sided P-value). Names are
#'       `omic1` ... `omicK`. Ready to pass directly to [ORBIT_P()],
#'       [ORBIT_Rank()], or [ORBIT_cor()].}
#'     \item{`Z_mat`}{The underlying `n` by `K` normal draws, with
#'       row names `feature1` ... `featureN` and column names
#'       `omic1` ... `omicK`.}
#'     \item{`is_signal`}{Logical vector of length `n`, named by
#'       feature. `TRUE` for up- or down-regulated features, `FALSE`
#'       for null features.}
#'     \item{`group`}{Character vector of length `n`, named by
#'       feature, with values `"null"`, `"up"`, or `"down"`
#'       indicating each feature's true class.}
#'     \item{`Sigma_null`}{The `K` by `K` exchangeable correlation
#'       matrix used for null features.}
#'     \item{`Sigma_signal`}{The `K` by `K` exchangeable correlation
#'       matrix used for signal features.}
#'   }
#'
#' @export
#'
#' @examples
#' # 80% null, 10% up, 10% down, with between-omic correlation
#' sim <- ORBIT_simulate_signal(
#'   n = 500, K = 3,
#'   pi_null = 0.8, pi_up = 0.1, pi_down = 0.1,
#'   mu_signal = 2,
#'   rho_null = 0.2, rho_signal = 0.5,
#'   seed = 1
#' )
#'
#' # Returned structure: list of per-omic data frames + ground-truth labels
#' names(sim)
#' head(sim$omics_list[[1]])
#' table(sim$is_signal)
ORBIT_simulate_signal <- function(n = 1000,
                                  K = 2,
                                  pi_null = 0.8,
                                  pi_up = 0.1,
                                  pi_down = 0.1,
                                  mu_signal = 2,
                                  rho_null = 0,
                                  rho_signal = 0,
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

  pis <- c(pi_null, pi_up, pi_down)
  if (!is.numeric(pis) || length(pis) != 3 ||
      any(!is.finite(pis)) || any(pis < 0) || any(pis > 1)) {
    stop("`pi_null`, `pi_up`, and `pi_down` must each lie in [0, 1].")
  }
  if (abs(sum(pis) - 1) > 1e-8) {
    stop("`pi_null + pi_up + pi_down` must equal 1.")
  }

  if (!is.numeric(mu_signal) || length(mu_signal) != 1 ||
      !is.finite(mu_signal)) {
    stop("`mu_signal` must be a single finite number.")
  }

  rho_lower <- -1 / (K - 1)
  for (nm in c("rho_null", "rho_signal")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
      stop(sprintf("`%s` must be a single finite number.", nm))
    }
    if (val <= rho_lower || val >= 1) {
      stop(sprintf(
        "`%s` must lie in the open interval (%.4f, 1) for K = %d.",
        nm, rho_lower, K
      ))
    }
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("`seed` must be a single number or NULL.")
    }
    set.seed(seed)
  }

  n <- as.integer(n)
  K <- as.integer(K)

  # ----------------------------------------------------------
  # 2. Build exchangeable correlation matrices and Cholesky factors
  # ----------------------------------------------------------
  build_exch <- function(rho, K) {
    S <- matrix(rho, nrow = K, ncol = K)
    diag(S) <- 1
    S
  }
  Sigma_null   <- build_exch(rho_null,   K)
  Sigma_signal <- build_exch(rho_signal, K)
  L_null   <- chol(Sigma_null)
  L_signal <- chol(Sigma_signal)

  # ----------------------------------------------------------
  # 3. Assign group labels and pre-allocate Z
  # ----------------------------------------------------------
  group <- sample(
    c("null", "up", "down"),
    size = n,
    replace = TRUE,
    prob = pis
  )
  is_signal <- group != "null"

  Z <- matrix(NA_real_, nrow = n, ncol = K)

  idx_null <- which(group == "null")
  idx_up   <- which(group == "up")
  idx_down <- which(group == "down")

  # ----------------------------------------------------------
  # 4. Draw Z for each component
  # ----------------------------------------------------------
  if (length(idx_null) > 0) {
    Z_iid <- matrix(stats::rnorm(length(idx_null) * K),
                    nrow = length(idx_null), ncol = K)
    Z[idx_null, ] <- Z_iid %*% L_null
  }

  n_signal <- length(idx_up) + length(idx_down)
  if (n_signal > 0) {
    # Draw correlated standard normals, then add the signed mean
    Z_iid_sig <- matrix(stats::rnorm(n_signal * K),
                        nrow = n_signal, ncol = K)
    Z_sig_zero_mean <- Z_iid_sig %*% L_signal

    if (length(idx_up) > 0) {
      Z[idx_up, ] <- Z_sig_zero_mean[seq_along(idx_up), , drop = FALSE] +
        mu_signal
    }
    if (length(idx_down) > 0) {
      down_rows <- seq.int(length(idx_up) + 1L, n_signal)
      Z[idx_down, ] <- Z_sig_zero_mean[down_rows, , drop = FALSE] -
        mu_signal
    }
  }

  # ----------------------------------------------------------
  # 5. Two-sided P-values and signs
  # ----------------------------------------------------------
  P <- 2 * stats::pnorm(-abs(Z))
  P <- pmax(P, .Machine$double.xmin)
  sign_mat <- ifelse(Z >= 0, 1, -1)

  # ----------------------------------------------------------
  # 6. Names
  # ----------------------------------------------------------
  feat_names <- paste0("feature", seq_len(n))
  omic_names <- paste0("omic",    seq_len(K))
  dn <- list(feat_names, omic_names)

  dimnames(P)            <- dn
  dimnames(sign_mat)     <- dn
  dimnames(Z)            <- dn
  dimnames(Sigma_null)   <- list(omic_names, omic_names)
  dimnames(Sigma_signal) <- list(omic_names, omic_names)
  names(group)           <- feat_names
  names(is_signal)       <- feat_names

  # ----------------------------------------------------------
  # 7. Build per-omic data frames
  # ----------------------------------------------------------
  omics_list <- setNames(lapply(seq_len(K), function(j) {
    data.frame(
      feature = feat_names,
      sign    = sign_mat[, j],
      stat    = P[, j],
      stringsAsFactors = FALSE
    )
  }), omic_names)

  list(
    omics_list   = omics_list,
    Z_mat        = Z,
    is_signal    = is_signal,
    group        = group,
    Sigma_null   = Sigma_null,
    Sigma_signal = Sigma_signal
  )
}
