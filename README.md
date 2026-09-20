
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ORBIT <img src="https://github.com/user-attachments/assets/15c5d57c-1d4e-48a0-a2da-621cd0bf181c" align="right" height="139" alt="" />

<!-- badges: start -->

<!-- badges: end -->

**ORBIT** integrates per-feature summary statistics across multiple omic
layers into a single calibrated P-value, while accounting for
between-omic correlation. It needs only the output of your per-omic
differential analyses — effect direction plus P-value or rank statistic
— not the raw sample-level data.

## Interactive website and examples

**[Click
here](https://yang-luo-lab.github.io/Rank-based-integration-identifies-convergent-disease-mechanisms-across-omics/)**
to explore the CKD and DCM results interactively, or run ORBIT on your
own summary statistics directly in the browser.

## Description

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

## System requirements

**Software dependencies**

- R (tested on version 4.5.2)
- R packages: `statmod` (tested on version 1.5.2) and `stats` (included
  with base R)
- Optional, for running tests: `testthat` (\>= 3.0.0)

**Operating systems**

ORBIT is written in pure R with no compiled code and is expected to run
on any platform supported by R (macOS, Windows, Linux).

**Versions tested**

- macOS Tahoe 26.6.1, R 4.5.2, statmod 1.5.2

**Hardware**

No non-standard hardware is required. ORBIT runs on a standard desktop
or laptop computer.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("yang-luo-lab/ORBIT")
```

Typical install time on a standard desktop computer is about 10 seconds.

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

**Expected output.** `res` is a data frame with one row per feature:
`Feature` (feature ID), `N` (number of omics in which the feature was
observed), `P` (combined P-value) and `Direction` (direction of the
combined effect, 1 or -1). With `seed = 1` the estimated `rho` is 0.3005
and the final confusion table matches the one shown above.

**Expected run time.** The full example above runs in under 1 second on
a standard desktop computer (0.1 s on the tested system).

## Instructions for use

**Running ORBIT on your own data**

ORBIT takes a named list of data frames, one per omic layer. Each data
frame has one row per feature and the following columns:

| Column | Type | Description |
|----|----|----|
| `feature` | character | Feature ID shared across omics (e.g. gene symbol) |
| `sign` | numeric | Direction of effect from the per-omic analysis, `1` or `-1` |
| `stat` | numeric | Per-feature P-value (for `ORBIT_P()`) or rank statistic (for `ORBIT_Rank()`) |
| `significance` | logical | Optional. Per-omic significance call (e.g. FDR \< 0.05), used by `ORBIT_cor()` with `sig_mode = "column"` to exclude signal features when estimating `rho` |

Features do not need to be present in every omic; ORBIT takes the union
of features and reports the number of contributing omics in `N`.

``` r
library(ORBIT)

omics_list <- list(
  rna     = data.frame(feature = rna_de$gene,  sign = sign(rna_de$logFC),  stat = rna_de$pvalue),
  protein = data.frame(feature = prot_de$gene, sign = sign(prot_de$logFC), stat = prot_de$pvalue)
)
direction <- c(rna = 1, protein = 1)

# Flag per-omic significant features so they are excluded from the rho estimate
omics_list$rna$significance     <- p.adjust(omics_list$rna$stat, "BH") < 0.05
omics_list$protein$significance <- p.adjust(omics_list$protein$stat, "BH") < 0.05

fit <- ORBIT_cor(omics_list, direction, sig_mode = "column")
res <- ORBIT_P(omics_list, direction, rho = fit$rho_mat)
res$FDR <- p.adjust(res$P, method = "BH")
```

ORBIT can also be run on your own summary statistics without installing
R, using the [interactive
website](https://yang-luo-lab.github.io/Rank-based-integration-identifies-convergent-disease-mechanisms-across-omics/).

**Reproducing the results in the manuscript**

Code and summary statistics for reproducing the CKD and DCM analyses are
available at
<https://github.com/yang-luo-lab/Rank-based-integration-identifies-convergent-disease-mechanisms-across-omics>.

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

## Citation

If you use ORBIT in your research, please cite:

> \[Author(s) (Year). Title. Journal, vol(num), pages.\]
