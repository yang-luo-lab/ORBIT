
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ORBIT

<!-- badges: start -->

<!-- badges: end -->

**ORBIT** integrates per-feature summary statistics across multiple omic
layers into a single calibrated P-value, while accounting for
between-omic correlation. It needs only the output of your per-omic
differential analyses — effect direction plus P-value or rank statistic
— not the raw sample-level data.

The package centers on three functions:

- **`ORBIT_cor()`** — estimate the between-omic correlation `ρ` from
  per-pair non-significant features (background-filtered,
  signal-resistant).
- **`ORBIT_Rank()`** — combine per-feature rank statistics across omics
  under a Variance-Gamma null parameterized by `ρ`.
- **`ORBIT_P()`** — combine per-feature P-values across omics, with
  optional auto-estimation of `ρ`.

`ρ` can be a scalar (uniform across pairs), a named symmetric matrix
(per-pair, e.g. from `ORBIT_cor()$rho_mat`), or left as `NULL` in
`ORBIT_P()` for automatic estimation.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zifengqiu/ORBIT")
```

## Quick example

Simulate three omics with 20% true signal, estimate the null
between-omic correlation, and combine evidence end-to-end:

``` r
library(ORBIT)

# Simulate 500 features × 3 omics; 80% null, 20% signal,
# rho_null = 0.3 in null features, rho_signal = 0.5 in signal features.
sim <- ORBIT_simulate_signal(
  n = 500, K = 3,
  pi_null = 0.8, pi_up = 0.1, pi_down = 0.1,
  mu_signal = 2,
  rho_null = 0.3, rho_signal = 0.5,
  seed = 1
)
direction <- setNames(rep(1, 3), names(sim$omics_list))

# Step 1. Estimate rho from background, filtering out signal features.
# In practice the `significance` flag is your per-omic DE call
# (e.g. FDR < 0.05). Here we use the simulated truth.
omics <- sim$omics_list
for (j in seq_along(omics))
  omics[[j]]$significance <- sim$is_signal[omics[[j]]$feature]

fit <- ORBIT_cor(omics, direction, sig_mode = "column", seed = 1)
#> Per-pair bg size: median=405, range=[405, 405] over 3 pairs
#> rho = 0.3005
fit$rho        # recovers rho_null ~ 0.3, unaffected by rho_signal
#> [1] 0.3004774

# Step 2. Combine evidence using the estimated per-pair rho matrix.
res <- ORBIT_P(sim$omics_list, direction, rho = fit$rho_mat)
#> Union of features across 3 omics: 500
#> Shared across ALL omics: 500 (100.0% of union)
#> Per-feature effective rho over 500 features with N>=2: median=0.3005, range=[0.3005, 0.3005]

# Step 3. FDR control + sanity check against the simulated truth.
fdr   <- p.adjust(res$P, method = "BH")
truth <- sim$is_signal[res$Feature]
table(predicted = fdr < 0.05, truth = truth)
#>          truth
#> predicted FALSE TRUE
#>     FALSE   401   65
#>     TRUE      4   30
```

## How it works

ORBIT models the combined statistic `D = Σ signed_score_j` under a
symmetric Variance-Gamma null parameterized by effective between-omic
correlation `ρ_eff`. The null is exact (not approximated as chi-squared,
unlike Brown’s method), evaluated via 32-node Gauss-Laguerre quadrature
with closed forms at `r_eff = 1` (Laplace) and `r_eff = 2`.

`ORBIT_cor()`’s background filter — excluding features significant in
either omic of a pair before estimating `ρ` — is the signal-resistance
step. Without it, concordant signals across omics push features toward
extreme ranks together and artificially inflate `ρ`. With it, the
estimate stays close to the null correlation regardless of signal
strength.

For full methodology, comparisons, and worked examples, see the
introductory vignette:

``` r
vignette("orbit-intro", package = "ORBIT")
```

## Citation

If you use ORBIT in your research, please cite:

> \[Author(s) (Year). Title. Journal, vol(num), pages.\]
