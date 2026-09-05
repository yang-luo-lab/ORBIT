#' Estimate Average Between-Omic Correlation from Per-Pair Non-Significant Features
#'
#' For each pair of omics, restricts to features that are present in
#' both omics and flagged as non-significant in both, recomputes ranks
#' within that subset, and forms a signed score whose Pearson
#' correlation across omics estimates the residual between-omic
#' correlation under the null. Returns the mean off-diagonal
#' correlation, the full pairwise correlation matrix, and a per-pair
#' diagnostic table.
#'
#' @param omics_list Named list of data frames, one per omic. Each must
#'   contain columns named by `feature_col`, `sign_col`, `stat_col`, and
#'   (when `sig_mode = "column"`) `sig_col`.
#' @param direction Named numeric vector of +1/-1; names must match
#'   `names(omics_list)`.
#' @param sig_mode `"column"` reads significance from `sig_col` (must
#'   be 0/1, no NA). `"thr"` ignores `sig_col` and flags the
#'   most-significant `floor(m * sig_thr)` features (minimum 1) per
#'   omic, ranking by `stat` with random tie-breaking controlled by
#'   `seed`. The direction of "most significant" is set by
#'   `stat_direction`.
#' @param sig_thr Required when `sig_mode = "thr"`; a number in `(0, 1)`.
#' @param stat_direction Convention for `stat_col`.
#'   `"smaller_better"` (default) treats smaller `stat` as more
#'   significant (p-value style); `"larger_better"` treats larger
#'   `stat` as more significant (`|t-stat|`, `-log10(p)`, enrichment
#'   score, etc.). Either way `stat` must be non-negative; the sign of
#'   the effect belongs in `sign_col`.
#' @param feature_col,sign_col,stat_col,sig_col Column names. Defaults
#'   `"feature"`, `"sign"`, `"stat"`, `"significance"`. `sig_col` is
#'   consulted only when `sig_mode = "column"`.
#' @param min_bg Minimum per-pair background size required to estimate
#'   that pair. Pairs below are left as `NA` in `rho_mat`. Default `30`.
#' @param bg_frac_tol Threshold on `bg_frac = n_bg / n_overlap` used
#'   for the low-bg-fraction warning. A pair is flagged when its
#'   `bg_frac < bg_frac_tol`. Default `0.5`. Set to `0` to silence.
#' @param seed Optional seed for random tie-breaking. The global RNG
#'   state is restored on exit.
#'
#' @details
#' **Score construction.** Within each pair's background, the score
#' for a feature in omic `k` is `sign * direction[k] * (-log(adj_rank))`,
#' where `adj_rank = r / (m_bg + 1)` and `r` is the rank of `stat`
#' arranged so that the most-significant feature has rank 1 (controlled
#' by `stat_direction`). Under the null `|score|` is approximately
#' `Exp(1)`-distributed and the signed score is marginally
#' Laplace(0, 1).
#'
#' **Per-pair background.** For pair `(i, j)` the background is the
#' set of features that are (a) observed in both omic `i` and omic `j`,
#' and (b) have `sig = 0` in both omic `i` and omic `j`. There is no
#' cross-pair "global null" set: features that are significant in some
#' third omic `k` are unaffected as long as they are non-significant in
#' both `i` and `j`. Ranks `r` are recomputed *within* each pair's
#' background before forming the score, so removing significant
#' features does not bias the rank distribution.
#'
#' **Aggregate `rho`.** `rho` is the arithmetic mean of the upper-
#' triangular finite pair correlations (no trimming applied).
#'
#' **Diagnostic table.** For each unordered pair, `diag` records:
#' overlap size (`n_overlap`), background size (`n_bg`), background
#' fraction (`bg_frac = n_bg / n_overlap`), Pearson correlation
#' (`rho`), a Tukey-IQR signed deviation (`tukey_dev`: distance to the
#' nearest IQR fence in IQR units, `0` if inside the box, sign
#' indicates direction; `|tukey_dev| > 1.5` is the classical Tukey
#' outlier threshold but no automatic exclusion is applied), and a
#' leave-one-out `rho_loo` (mean of the other valid pair correlations,
#' or the overall mean of valid pairs for rows whose own correlation
#' is `NA`). Rows are sorted by `|tukey_dev|` descending so the most
#' unusual pairs appear first. `tukey_dev` is `NA` when fewer than
#' four valid pairs exist or the IQR is zero; `rho_loo` is `NA` when
#' fewer than two valid pairs exist.
#'
#' **Messages and warnings.** The function emits two `message()`s
#' giving the per-pair background size summary (median, range, count)
#' and the value of `rho`. A single `warning()` is issued when any
#' pair has `n_bg < min_bg` or `bg_frac < bg_frac_tol`; specific pairs
#' are identified in `$diag`.
#'
#' @return A list with:
#' \itemize{
#'   \item `rho`: mean of the upper-triangle of `rho_mat` over valid
#'     pairs (no trimming). `NA` if no pair has a finite correlation.
#'   \item `rho_mat`: pairwise Pearson correlation on per-pair
#'     background. Diagonal is `1`; entries are `NA` when the per-pair
#'     background is below `min_bg` or the correlation is undefined.
#'   \item `diag`: data frame with one row per unordered pair and
#'     columns `omic_i`, `omic_j`, `n_overlap`, `n_bg`, `bg_frac`,
#'     `rho`, `tukey_dev`, `rho_loo`. Sorted by `|tukey_dev|`
#'     descending.
#'   \item `n_total`: size of the feature union across omics.
#'   \item `tukey_quartiles`: named vector `c(Q1, Q3, IQR)` of the
#'     pair-correlation distribution; entries are `NA` when fewer than
#'     four valid pairs exist.
#'   \item `settings`: list of resolved values `min_bg` and
#'     `bg_frac_tol`.
#' }
#' @importFrom stats cor median quantile setNames
#' @examples
#' # Simulate data with real signal: 80% null, 20% true signal,
#' # rho_null = 0.3 (null between-omic correlation we want to recover),
#' # rho_signal = 0.5 (extra correlation among signal features).
#' sim <- ORBIT_simulate_signal(
#'   n = 500, K = 3,
#'   pi_null = 0.8, pi_up = 0.1, pi_down = 0.1,
#'   mu_signal = 2,
#'   rho_null = 0.3, rho_signal = 0.5,
#'   seed = 1
#' )
#'
#' # Attach a `significance` flag per omic. In practice this is your
#' # DE call (e.g. FDR < 0.05); here we use the simulated truth.
#' omics <- sim$omics_list
#' for (j in seq_along(omics))
#'   omics[[j]]$significance <- sim$is_signal[omics[[j]]$feature]
#' direction <- setNames(rep(1, 3), names(omics))
#'
#' fit <- ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
#' fit$rho           # recovers rho_null ~ 0.3 (signal filtered out)
#' fit$rho_mat       # pairwise matrix; feed into ORBIT_Rank / ORBIT_P
#' fit$diag          # per-pair diagnostics
#' @export
ORBIT_cor <- function(omics_list,
                      direction,
                      sig_mode       = c("column", "thr"),
                      sig_thr        = NULL,
                      stat_direction = c("smaller_better", "larger_better"),
                      feature_col    = "feature",
                      sign_col       = "sign",
                      stat_col       = "stat",
                      sig_col        = "significance",
                      min_bg         = 30L,
                      bg_frac_tol    = 0.5,
                      seed           = NULL) {

  # ---- 0. Match args + parameter checks ----
  sig_mode       <- match.arg(sig_mode)
  stat_direction <- match.arg(stat_direction)

  check_col <- function(x, nm) {
    if (!is.character(x) || length(x) != 1 || is.na(x) || x == "")
      stop("`", nm, "` must be a single non-empty character string.")
  }
  check_col(feature_col, "feature_col"); check_col(sign_col, "sign_col")
  check_col(stat_col,    "stat_col")
  if (sig_mode == "column") check_col(sig_col, "sig_col")

  cols_used <- c(feature_col, sign_col, stat_col,
                 if (sig_mode == "column") sig_col)
  if (anyDuplicated(cols_used))
    stop("Column-name parameters must all be different.")

  if (sig_mode == "thr") {
    if (is.null(sig_thr) || !is.numeric(sig_thr) || length(sig_thr) != 1 ||
        !is.finite(sig_thr) || sig_thr <= 0 || sig_thr >= 1)
      stop("`sig_thr` must be a single number in (0, 1) when sig_mode = 'thr'.")
  }
  if (!is.numeric(min_bg) || length(min_bg) != 1 || min_bg < 2)
    stop("`min_bg` must be an integer >= 2.")
  min_bg <- as.integer(min_bg)

  if (!is.numeric(bg_frac_tol) || length(bg_frac_tol) != 1 ||
      !is.finite(bg_frac_tol) || bg_frac_tol < 0 || bg_frac_tol > 1)
    stop("`bg_frac_tol` must be a single number in [0, 1].")

  sig_rank <- function(t) {
    if (stat_direction == "larger_better") rank(-t, ties.method = "random")
    else                                    rank( t, ties.method = "random")
  }

  # ---- 1. omics_list structural checks ----
  if (!is.list(omics_list) || is.null(names(omics_list)) ||
      any(names(omics_list) == "") || anyDuplicated(names(omics_list)))
    stop("`omics_list` must be a named list with unique, non-empty omic names.")
  omic_names <- names(omics_list); n_omic <- length(omics_list)
  if (n_omic < 2L) stop("`omics_list` must contain at least 2 omics.")

  # ---- 2. Per-omic checks; standardize ----
  for (j in seq_len(n_omic)) {
    df_j <- omics_list[[j]]; nm <- omic_names[j]
    if (!is.data.frame(df_j)) stop("Omic '", nm, "' is not a data frame.")
    miss <- setdiff(cols_used, colnames(df_j))
    if (length(miss))
      stop("Omic '", nm, "' missing column(s): ", paste(miss, collapse = ", "))
    f <- df_j[[feature_col]]; s <- df_j[[sign_col]]; t <- df_j[[stat_col]]
    if (is.factor(f)) f <- as.character(f)
    if (!is.character(f))
      stop("Omic '", nm, "': `", feature_col, "` must be character/factor.")
    if (anyNA(f) || any(f == ""))
      stop("Omic '", nm, "': `", feature_col, "` has missing/empty values.")
    if (anyDuplicated(f))
      stop("Omic '", nm, "': duplicate values in `", feature_col, "`.")
    if (!is.numeric(s)) stop("Omic '", nm, "': `", sign_col, "` not numeric.")
    if (anyNA(s))      stop("Omic '", nm, "': `", sign_col, "` has NA.")
    if (any(!s %in% c(-1, 1)))
      stop("Omic '", nm, "': `", sign_col, "` must be -1 or +1.")
    if (!is.numeric(t)) stop("Omic '", nm, "': `", stat_col, "` not numeric.")
    if (anyNA(t))       stop("Omic '", nm, "': `", stat_col, "` has NA.")
    if (any(t < 0))
      stop("Omic '", nm, "': `", stat_col, "` has negative values.")
    if (sig_mode == "column") {
      sig <- df_j[[sig_col]]
      if (!is.numeric(sig) && !is.logical(sig))
        stop("Omic '", nm, "': `", sig_col, "` must be numeric/logical (0/1).")
      sig <- as.integer(sig)
      if (anyNA(sig)) stop("Omic '", nm, "': `", sig_col, "` has NA.")
      if (!all(sig %in% c(0L, 1L)))
        stop("Omic '", nm, "': `", sig_col, "` must contain only 0/1.")
    } else sig <- NA_integer_
    omics_list[[j]] <- data.frame(feature = f, sign = s, stat = t,
                                  significance = sig,
                                  stringsAsFactors = FALSE)
  }

  # ---- 3. direction ----
  if (!is.numeric(direction) || is.null(names(direction)))
    stop("`direction` must be a named numeric vector.")
  if (!setequal(names(direction), omic_names))
    stop("Names of `direction` must match names of `omics_list`.")
  direction <- direction[omic_names]
  if (anyNA(direction) || !all(direction %in% c(-1, 1)))
    stop("`direction` must contain only -1 or +1.")

  # ---- 4. Seed handling ----
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

  # ---- 5. thr mode: per-omic top-rank significance ----
  if (sig_mode == "thr") {
    for (j in seq_len(n_omic)) {
      t   <- omics_list[[j]]$stat; m_j <- length(t)
      if (m_j == 0) next
      n_top <- max(1L, floor(m_j * sig_thr))
      omics_list[[j]]$significance <- as.integer(sig_rank(t) <= n_top)
    }
  }

  # ---- 6. Feature universe ----
  feature_lists <- lapply(omics_list, `[[`, "feature")
  all_features  <- unique(unlist(feature_lists, use.names = FALSE))
  n_total       <- length(all_features)

  # ---- 7. Wide matrices on the union ----
  has_feat  <- matrix(FALSE,    n_total, n_omic,
                      dimnames = list(all_features, omic_names))
  sig_wide  <- matrix(0L,       n_total, n_omic,
                      dimnames = list(all_features, omic_names))
  stat_wide <- matrix(NA_real_, n_total, n_omic,
                      dimnames = list(all_features, omic_names))
  sign_wide <- matrix(NA_real_, n_total, n_omic,
                      dimnames = list(all_features, omic_names))
  for (j in seq_len(n_omic)) {
    df_j <- omics_list[[j]]; if (!nrow(df_j)) next
    idx  <- match(df_j$feature, all_features)
    has_feat [idx, j] <- TRUE
    sig_wide [idx, j] <- df_j$significance
    stat_wide[idx, j] <- df_j$stat
    sign_wide[idx, j] <- df_j$sign
  }

  # ---- 8. Pairwise: bg = overlap & non-sig in BOTH omics; rerank; Pearson ----
  rho_mat        <- matrix(NA_real_, n_omic, n_omic,
                           dimnames = list(omic_names, omic_names))
  n_pair         <- matrix(0L, n_omic, n_omic,
                           dimnames = list(omic_names, omic_names))
  n_pair_overlap <- matrix(0L, n_omic, n_omic,
                           dimnames = list(omic_names, omic_names))
  diag(rho_mat) <- 1

  for (i in seq_len(n_omic - 1)) for (j in seq.int(i + 1, n_omic)) {
    overlap_ij <- has_feat[, i] & has_feat[, j]
    n_pair_overlap[i, j] <- n_pair_overlap[j, i] <- sum(overlap_ij)
    bg_ij   <- overlap_ij & sig_wide[, i] == 0L & sig_wide[, j] == 0L
    n_bg_ij <- sum(bg_ij)
    n_pair[i, j] <- n_pair[j, i] <- n_bg_ij

    if (n_bg_ij < min_bg) next

    stat_i <- stat_wide[bg_ij, i]; sign_i <- sign_wide[bg_ij, i]
    stat_j <- stat_wide[bg_ij, j]; sign_j <- sign_wide[bg_ij, j]
    m_ij   <- n_bg_ij
    adj_i  <- sig_rank(stat_i) / (m_ij + 1)
    adj_j  <- sig_rank(stat_j) / (m_ij + 1)
    score_i <- sign_i * direction[i] * (-log(adj_i))
    score_j <- sign_j * direction[j] * (-log(adj_j))
    r_t <- suppressWarnings(stats::cor(score_i, score_j, method = "pearson"))
    if (is.finite(r_t)) rho_mat[i, j] <- rho_mat[j, i] <- r_t
  }

  # ---- 9. Per-pair diagnostic table + summary ----
  ut       <- upper.tri(rho_mat)
  pair_idx <- which(ut, arr.ind = TRUE)
  cor_vec  <- rho_mat[ut]
  finite_i <- is.finite(cor_vec)
  n_valid  <- sum(finite_i)
  tukey_dev <- rep(NA_real_, length(cor_vec))
  q1 <- q3 <- iq <- NA_real_
  if (n_valid >= 4) {
    qq <- stats::quantile(cor_vec[finite_i], c(0.25, 0.75), names = FALSE)
    q1 <- qq[1]; q3 <- qq[2]; iq <- q3 - q1
    if (iq > 0) {
      above <- pmax(0, cor_vec - q3) / iq
      below <- pmax(0, q1 - cor_vec) / iq
      tukey_dev <- above - below
    }
  }
  rho <- if (n_valid) mean(cor_vec[finite_i]) else NA_real_
  if (!is.finite(rho)) rho <- NA_real_
  rho_loo <- rep(NA_real_, length(cor_vec))
  if (n_valid > 1) {
    cor_sum_valid      <- sum(cor_vec[finite_i])
    rho_loo[finite_i]  <- (cor_sum_valid - cor_vec[finite_i]) / (n_valid - 1)
    rho_loo[!finite_i] <- cor_sum_valid / n_valid
  }

  diag_df <- data.frame(
    omic_i    = omic_names[pair_idx[, 1]],
    omic_j    = omic_names[pair_idx[, 2]],
    n_overlap = as.integer(n_pair_overlap[ut]),
    n_bg      = as.integer(n_pair[ut]),
    bg_frac   = n_pair[ut] / pmax(n_pair_overlap[ut], 1),
    rho       = cor_vec,
    tukey_dev = tukey_dev,
    rho_loo   = rho_loo,
    stringsAsFactors = FALSE
  )
  ord <- order(-abs(replace(diag_df$tukey_dev, is.na(diag_df$tukey_dev), 0)))
  diag_df <- diag_df[ord, , drop = FALSE]
  rownames(diag_df) <- NULL

  bg_vec     <- n_pair[ut]
  n_small_bg <- sum(bg_vec < min_bg)
  n_low_frac <- sum(diag_df$bg_frac < bg_frac_tol, na.rm = TRUE)

  message(sprintf("Per-pair bg size: median=%d, range=[%d, %d] over %d pairs",
                  as.integer(stats::median(bg_vec)),
                  min(bg_vec), max(bg_vec), length(bg_vec)))
  if (is.finite(rho))
    message(sprintf("rho = %.4f", rho))

  if (n_small_bg > 0 || n_low_frac > 0)
    warning(sprintf(
      "%d pair(s) with bg < %d; %d pair(s) with bg fraction < %.0f%% of overlap. Inspect `$diag` to identify them.",
      n_small_bg, min_bg, n_low_frac, 100 * bg_frac_tol))

  # ---- 10. Return ----
  list(
    rho             = rho,
    rho_mat         = rho_mat,
    diag            = diag_df,
    n_total         = n_total,
    tukey_quartiles = c(Q1 = q1, Q3 = q3, IQR = iq),
    settings        = list(min_bg      = min_bg,
                           bg_frac_tol = bg_frac_tol)
  )
}


