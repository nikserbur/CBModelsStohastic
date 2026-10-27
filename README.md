# CBModelsStohastic

[Русская версия](README.ru.md)

Research code in R and C++ for fat-tailed stochastic extensions of interest-rate modelling for Russia: tail diagnostics, a fat-tail adjustment of R², HAR volatility, hidden Markov regimes and an augmented Taylor rule with Monte Carlo scenarios for the RUONIA rate.

Companion code for the article [*Stochastic modelling of monetary policy: a modular architecture for analysing extreme scenarios*](http://scipro.ru/article/07-03-2026) (Burakov, International Journal of Professional Science, 2026, No. 3-1, pp. 47–66).

## What is in this repository

The article describes a modular pipeline:
- latent factors from copula normalisation and SVD;
- HAR / stochastic volatility with calibrated tail distributions;
- a hidden Markov model (HMM) for regime switching;
- an augmented Taylor rule as the policy reaction operator;
- Monte Carlo fan charts with tail-risk metrics (VaR / CTE);
- evaluation with Taleb's R².

This repository contains only part of that pipeline:

| Component | Code in this repository | Status |
|---|---|---|
| Tail diagnostics of distributions | [`1.R/3.dist_bust.R`](1.R/3.dist_bust.R), [`1.R/10.DashboardData.R`](1.R/10.DashboardData.R) | Runs standalone on simulated data |
| Fat-tail adjustment of R² (Taleb) | [`1.R/4.LmR2BoostingStohastic.R`](1.R/4.LmR2BoostingStohastic.R) | Runs standalone on simulated data |
| HAR volatility, volatility of volatility, GARCH(1,1) benchmark | [`1.R/Demo_Volatility_Analysis.R`](1.R/Demo_Volatility_Analysis.R) | Runs standalone on simulated data |
| RDC-based feature selection (C++) | [`2.Cpp/RDCScpp.cpp`](2.Cpp/RDCScpp.cpp) | Compiles with Rcpp; no R script calls it yet |
| HMM regimes (5 states) | [`1.R/5.HiddenMarkingsModels.R`](1.R/5.HiddenMarkingsModels.R) | Needs functions that are not in the repository |
| SVD forecasting | [`1.R/Demo_SVD_Forecasting.R`](1.R/Demo_SVD_Forecasting.R) | Needs `PredictTimeSeriesWithSVD()`, which is not in the repository |
| Augmented Taylor rule, tail calibration, Monte Carlo, VaR / CVaR | [`1.R/Main.R`](1.R/Main.R) | Draft: runs through data preparation, then stops |
| Copula normalisation | — | Not in the repository |

None of the scripts that use the real data in [`3.Data/`](3.Data) currently runs to completion. The empirical results in the article therefore cannot be reproduced from this repository alone. Details are in [Run order](#run-order) and [Scripts](#scripts).

## How to cite

If you use this code, please cite the article:

> Burakov N. S. Stochastic modelling of monetary policy: a modular architecture for analysing extreme scenarios. *International Journal of Professional Science*, 2026, No. 3-1, pp. 47–66. EDN: TPCNIM. URL: http://scipro.ru/article/07-03-2026

```bibtex
@article{burakov2026stochastic,
  author  = {Burakov, Nikita S.},
  title   = {Stochastic modelling of monetary policy: a modular architecture for analysing extreme scenarios},
  journal = {International Journal of Professional Science},
  year    = {2026},
  number  = {3-1},
  pages   = {47--66},
  note    = {EDN: TPCNIM},
  url     = {http://scipro.ru/article/07-03-2026}
}
```

## Requirements

- **R ≥ 4.0.** The four standalone scripts and the C++ example below were tested with R 4.4.1 on macOS.
- **R packages.** [`1.R/9.libraries.r`](1.R/9.libraries.r) lists every package the scripts use. Sourcing it installs any missing package from CRAN and attaches all of them.
  - `4.LmR2BoostingStohastic.R` and `Demo_Volatility_Analysis.R` attach their own packages with `library()`.
  - `qs` is optional. It is used only by the `save_qs_path` argument of `BuildDistDashboard()`. CRAN archived it on 2026-01-17, so install it from the CRAN archive if you need that argument.
- **C++ toolchain** (only for [`2.Cpp/RDCScpp.cpp`](2.Cpp/RDCScpp.cpp)). `Rcpp::sourceCpp()` needs a C++ compiler:
  - macOS: Xcode Command Line Tools (`xcode-select --install`);
  - Windows: the Rtools version that matches your R version;
  - Linux: `g++` and the R development headers (for example, `build-essential` and `r-base-dev` on Debian/Ubuntu).
- **Working directory.** All file paths in the scripts are relative to the repository root. Open [`CB_models_stohastic.Rproj`](CB_models_stohastic.Rproj) in RStudio, or call `setwd()` with the repository root.

## Data

[`3.Data/`](3.Data) contains the input files. `Main.R` gives https://www.cbr.ru/dkp/statistics/ as the source of its inputs. The article (section "Materials and methods") lists these sources for its monthly panel:
- Bank of Russia: RUONIA and the key rate;
- Rosstat and Bank of Russia surveys: inflation and inflation expectations;
- Rosstat and the Ministry of Finance: GDP and the credit impulse;
- Moscow Exchange: market volatility.

The repository does not record the source of each file, except where a file names its publisher in its own header.

| File | Content | Frequency | Period | Used by |
|---|---|---|---|---|
| `RC_F11_01_2010_T04_07_2025.xlsx` | RUONIA: rate (`ruo`), volume (`vol`) and related fields | Daily | 11 Jan 2010 – 4 Jul 2025 | `Main.R` (monthly mean of `ruo`) |
| `Инфляция и ключевая ставка Банка России_F01_12_2013_T18_06_2025.xlsx` | Key rate (% p.a.), CPI inflation (% y/y), inflation target | Monthly | Dec 2013 – May 2025 | `Main.R`, `5.HiddenMarkingsModels.R` |
| `кредитный_импульс_помесячно.csv` | Credit impulse | Monthly | Feb 2010 – Dec 2024 | `Main.R`, `5.HiddenMarkingsModels.R` |
| `vvp_pomesyachno.csv` | GDP level (`ВВП_месяц`) | Monthly | Jan 2011 – Mar 2025 | `Main.R`, `5.HiddenMarkingsModels.R` |
| `ожидаемая_инфляция.csv` | Expected inflation, % | Irregular until Mar 2020, monthly afterwards | Dec 2009 – Jun 2025 | `Main.R`, `5.HiddenMarkingsModels.R` |
| `budget.xlsx` | Budget series (`Value`), called the "budget impulse" in `Main.R` | Monthly | Jan 2010 – Mar 2025 | `Main.R` |
| `assets.xlsx` | MOEX and RTS index levels | Daily | Jan 2017 – May 2025 | `Main.R` (30-day rolling volatility of MOEX + RTS) |
| `sunctions.csv` | Sanctions indicator: unnamed columns, non-decreasing values, called "Скор санкций" in `Main.R` | Irregular | 21 Feb 2022 – 15 Jan 2025 | Read by `Main.R`, but not merged (the merge line is commented out) |
| `stavka.xlsx` | RUONIA export with the same layout as `RC_F11_…` | Daily | 1 Nov 2010 – 16 Jun 2025 | Unused |
| `VVP_kvartal.xlsx` | GDP tables: current prices, constant prices, seasonally adjusted, volume indices, deflators | Quarterly | 1995 – 2025 | Unused |
| `add_info.xlsx` | Bank of Russia enterprise monitoring survey (30 sheets by sector) | Monthly | Jan 2019 – Apr 2025 | Unused |
| `bal_of_payments_ap.xlsx` | Balance of payments of the Russian Federation, analytical presentation (USD million) | Quarterly, annual | Q1 1994 – Q1 2025 | Unused |
| `bb.xlsx` | Bank of Russia Statistics Department: financial condition statistics of non-financial organisations | — | Header dated Dec 2024 | Unused |
| `exchange_rate.xlsx` | Main derived indicators of rouble exchange-rate dynamics | Monthly, quarterly | Jan 2005 – May 2025 | Unused |
| `inflation.xlsx` | CPI growth, % m/m, seasonally adjusted, by component | Monthly | From Jan 2002 | Unused |

The monthly panel built in `Main.R` covers Feb 2017 – May 2025 (100 months). It starts in 2017 because the merge is anchored on the MOEX/RTS volatility series.

**File names.** Three data files have Cyrillic names, and the scripts read them by literal path. On systems whose native encoding is not UTF-8 (for example, R versions before 4.2 on Windows), these paths may fail to resolve. The files are not renamed because the scripts refer to them by name.

## Run order

Run these commands in R from the repository root. `print.eval = TRUE` displays the plots that the scripts create at top level.

```r
# 0. Install (if needed) and attach the packages listed in 1.R/9.libraries.r
source("1.R/9.libraries.r")

# 1. Tail diagnostics on simulated data: normal vs stochastic volatility vs mixture
source("1.R/3.dist_bust.R", print.eval = TRUE)

# 2. Standard R² vs fat-tail-adjusted R² for regressions with fat-tailed residuals
source("1.R/4.LmR2BoostingStohastic.R", print.eval = TRUE)

# 3. Distribution dashboard: moments, normality tests, Hill and power-law tail estimates
source("1.R/10.DashboardData.R")
pareto_example$plot    # one of six examples; see also normal_example, cauchy_example, ...
pareto_example$stats

# 4. HAR volatility forecast with shocks, stochastic band and GARCH(1,1) benchmark (simulated data)
source("1.R/Demo_Volatility_Analysis.R")

# 5. Optional: compile the C++ RDC feature selection (see "RDC feature selection" below)
Rcpp::sourceCpp("2.Cpp/RDCScpp.cpp")
```

The following scripts **do not run to completion** in the current state of the repository:

- **`1.R/8.SourceAll.R`** sources `1.R/9.libraries.r` and four function files that are not in the repository: `HarVolatilityFunctions.R`, `HmmFunctions.R`, `BoostingFunctions.R` and `TalebR2Functions.R` in `1.R/1.functions/`. It reports each missing file and ends with a warning.
- **`1.R/Main.R`** loads and merges the data in `3.Data/` (section 1). It stops in section 2, because `stable_fit()` is not defined in the repository or in any package it loads. Later sections also depend on code that does not exist yet:
  - an undefined data object `dt`;
  - `estimate_expectations()`, from a `3_ml_expectations.R` that is not in the repository;
  - the columns `pi_target`, `volatility_index` and `sanctions_score`, which the data step does not create;
  - `pit_sample()`, which `scoringRules` does not provide.
- **`1.R/5.HiddenMarkingsModels.R`** loads the data, then stops at `InitializeHmmParameters()`. It also needs `EstimateHmmParameters()` and `ViterbiDecoding()`.
- **`1.R/Demo_SVD_Forecasting.R`** sources [`1.R/1.functions/all_fun.R`](1.R/1.functions/all_fun.R), which is empty, and calls `PredictTimeSeriesWithSVD()`, which is not defined. It also needs `data.table` and `plotly` to be attached first.

## Scripts

```
CBModelsStohastic/
├── CB_models_stohastic.Rproj        RStudio project; opening it sets the working directory to the repository root
├── README.md, README.ru.md          Documentation (English, Russian)
├── 1.R/
│   ├── 1.functions/all_fun.R        Placeholder for shared functions (empty)
│   ├── 3.dist_bust.R                Tail diagnostics: normal vs stochastic-volatility vs mixture samples
│   ├── 4.LmR2BoostingStohastic.R    Standard vs fat-tail-adjusted R² (Taleb) on simulated regressions
│   ├── 5.HiddenMarkingsModels.R     5-regime HMM for the key rate (needs missing functions)
│   ├── 8.SourceAll.R                Loads 9.libraries.r and the shared function files; reports missing ones
│   ├── 9.libraries.r                Installs and attaches all R packages
│   ├── 10.DashboardData.R           BuildDistDashboard(): moments, normality tests, Hill / power-law tails
│   ├── Demo_SVD_Forecasting.R       SVD forecast of simulated macro series (needs PredictTimeSeriesWithSVD())
│   ├── Demo_Volatility_Analysis.R   HAR volatility with shocks and vol-of-vol band; GARCH(1,1) benchmark
│   └── Main.R                       Draft: data panel, tails, augmented Taylor rule, Monte Carlo, VaR / CVaR
├── 2.Cpp/
│   └── RDCScpp.cpp                  cancor_cpp() and RDCSelection_cpp(): RDC-based feature selection (Rcpp)
└── 3.Data/                          Input data (see Data)
```

The file name `5.HiddenMarkingsModels.R` is kept as is: it means "Hidden Markov Models".

## Method

This section describes what the code implements. Symbols are used consistently throughout:

| Symbol | Meaning |
|---|---|
| $\nu$ | Tail index / Student-t degrees of freedom (`ALPHA_TAIL`, `alpha_tail` in the code) |
| $\alpha$ | Stability index of an α-stable law (`Main.R`, `10.DashboardData.R`) |
| $c$ | Confidence level of VaR / CVaR (`Main.R` uses $c = 0.95$) |
| $K$ | Number of HMM states ($K = 5$) |
| $P$ | HMM transition matrix |

### Tail diagnostics (`3.dist_bust.R`)

The script draws $N = 100\,000$ observations from three processes with base volatility $\sigma_0 = 0.3$:
- **normal:** $X \sim N(0, \sigma_0^2)$;
- **stochastic volatility:** $X \mid s \sim N(0, s^2)$, with $s = \sigma_0 e^{\eta}$ and $\eta \sim N(0, 0.6^2)$;
- **mixture:** with probability 0.01, $X \sim N(0, (5\sigma_0)^2)$; otherwise $X$ is a stochastic-volatility draw.

For each sample it reports:
- $P(|X| > k\sigma_0)$ for $k = 2, \dots, 5$; the thresholds are in units of $\sigma_0$, not of the sample standard deviation;
- the normal benchmark $2(1 - \Phi(k))$;
- the kurtosis (Pearson definition, equal to 3 for a normal distribution), mean and variance;
- QQ plots against the normal distribution and histograms of $\log|X|$.

### Fat-tail adjustment of R² (`4.LmR2BoostingStohastic.R`)

**Simulated model.** The script simulates $n = 10\,000$ observations of

$$y = 2 + 1.5\,x_1 - 0.8\,x_2 + \varepsilon, \qquad x_1 \sim N(0, 1), \quad x_2 \sim N(0, 1.5^2),$$

with three residual types:
- **normal:** $\varepsilon \sim N(0, 0.5^2)$;
- **fat-tailed:** $\varepsilon = 0.45\,t_\nu$ with $\nu = 2.8$, where 4 % of the draws are multiplied by 3.5;
- **mixture:** with probability 0.07, a $1.5\,t_3$ draw; otherwise $N(0, 0.5^2)$.

**R² variants.** After OLS, `CalculateRobustR2(y, ŷ, ν)` computes three R² values. Here $\hat\varepsilon_i = y_i - \hat y_i$:

$$R^2_{\text{std}} = 1 - \frac{\sum_i \hat\varepsilon_i^2}{\sum_i (y_i - \bar y)^2},$$

$$R^2_{\text{Taleb}} = \frac{\widehat{\operatorname{Var}}(\hat y)}{\widehat{\operatorname{Var}}(\hat y) + \frac{\nu}{\nu - 2}\,\widehat{\operatorname{Var}}(\hat\varepsilon)}, \qquad R^2_{\text{MLE}} = \max\left(0,\; 1 - \frac{\nu}{\nu - 2}\cdot\frac{\sum_i \hat\varepsilon_i^2}{\sum_i (y_i - \bar y)^2}\right).$$

- Both adjusted values are set to 0 when $\nu \le 2.1$, and $R^2_{\text{Taleb}}$ is clipped to $[0, 1]$.
- The factor $\nu/(\nu-2)$ is the variance of a standard Student-t variable with $\nu$ degrees of freedom.
- The script uses $\nu = 30$ for normal residuals, $\nu = 2.8$ for fat-tailed residuals and $\nu = \max(4,\, 0.8 \cdot 2.8) = 4$ for the mixture.

**Relation to Taleb (2020).** The code comments refer to Taleb's expression $R^2 = \dfrac{a^2(\nu-2)}{a^2(\nu-2) + \nu}$. It holds for $y = a x + \varepsilon$ with $\operatorname{Var}(x) = 1$ and a standard Student-t error $\varepsilon \sim t_\nu$, because then

$$R^2 = \frac{a^2}{a^2 + \nu/(\nu-2)}.$$

The code differs in one respect. It multiplies the *empirical* residual variance, which already includes the factor $\nu/(\nu-2)$, by $\nu/(\nu-2)$ again. Its adjusted values are therefore more conservative (lower) than the population formula. The script's plot title calls them "conservative R² correction methods".

**Hill estimate.** The script also reports the Hill estimate of the tail index of $|\hat\varepsilon|$ from the $k = \min(100, \lfloor n/4 \rfloor)$ largest values $X_{(1)} \ge X_{(2)} \ge \dots$:

$$\hat\nu_{\text{Hill}} = \left(\frac{1}{k} \sum_{i=1}^{k} \ln X_{(i)} - \ln X_{(k+1)}\right)^{-1}.$$

### Distribution dashboard (`10.DashboardData.R`)

`BuildDistDashboard(x)` returns a list with:
- a six-panel plotly figure (histogram with kernel and normal densities, ECDF, normal QQ plot, log-log survival plot, Hill plot, violin/box);
- a table of statistics: moments, quantiles, median, MAD, Jarque–Bera, Shapiro–Wilk (on at most 5,000 points) and Anderson–Darling p-values;
- right-tail estimates:
  - the median of the Hill estimates over a range of $k$ within the top `tail_fraction` of the data;
  - a power-law exponent from a regression of $\ln(1 - \hat F(x))$ on $\ln x$;
- a verbal tail label.

The script runs it on six simulated samples: normal, a mixture of normals, Pareto, normal + Pareto, Cauchy, and α-stable with $\alpha = 0.5$.

### HAR volatility (`Demo_Volatility_Analysis.R`)

**Data.** The script simulates 2,500 daily returns. Volatility follows a GARCH-type recursion with four crisis windows, calendar effects and a leverage effect. Innovations switch to Student-t when the volatility multiplier is high, and 15 flash-crash days are added. The first 2,200 days are used for training and the last 300 for testing.

**Realised volatility proxy:**

$$RV_t = \sqrt{252}\,\sum_{i=0}^{21} r_{t-i}^2.$$

**HAR on log-RV** (Corsi, 2009). `highfrequency::HARmodel()` fits the model with horizons $h \in \{1, 5, 22, 91, 365, 730\}$ days:

$$\ln RV_{t+1} = \beta_0 + \sum_{h} \beta_h\, \overline{\ln RV}{}_t^{(h)} + u_{t+1}, \qquad \overline{\ln RV}{}_t^{(h)} = \frac{1}{h}\sum_{i=0}^{h-1} \ln RV_{t-i}.$$

Corsi's original specification uses $h \in \{1, 5, 22\}$. The code comments attribute the log transformation to Clements and Preve.

**Forecast.** The forecast is built one step at a time over 300 days. At each step the model is re-estimated on the extended series and transformed back with $\exp(\cdot)$. The forecast is then multiplied by heuristic factors that are set in the code, not estimated:
- random shocks, whose probability combines the historical shock frequency, a 350-day cycle, calendar effects and the deviation of recent log-volatility;
- regime adjustments between 0.95 and 1.2.

**Stochastic band.** The script simulates 1,000 paths of a mean-reverting volatility process around the HAR forecast. The code calls it Heston-type and attributes it to Lu (2022), for which it gives no full reference. The parameters are $\kappa = 0.5$, $\sigma_v = 0.3$, leverage $-0.5$, and the weight of the HAR forecast is $0.6 + 0.2e^{-t/5}$. It reports quantiles from 1 % to 99 %.

**Benchmark and accuracy.** The benchmark is a GARCH(1,1) model with normal innovations (`rugarch`), annualised as $\hat\sigma_t\sqrt{252}$. The script compares MSE, MAE and MAPE on the test window.

### HMM regimes (`5.HiddenMarkingsModels.R`; needs missing functions)

**Data.** The monthly key rate $r_t$ from January 2014, together with:
- CPI inflation;
- the credit impulse;
- GDP;
- expected inflation.

Missing values are filled with `zoo::na.spline()`. The script derives:
- the inflation gap to the 4 % target;
- the key rate minus inflation;
- GDP growth y/y;
- the credit impulse relative to GDP;
- the monthly key-rate change.

**Model.** This is a Markov-switching model in the sense of Hamilton (1989). It has $K = 5$ hidden states $S_t$ and a transition matrix $P_{jk} = P(S_{t+1} = k \mid S_t = j)$. The script expects per-state emission means $\mu_k$ and standard deviations $\sigma_k$, that is, Gaussian emissions $r_t \mid S_t = k \sim N(\mu_k, \sigma_k^2)$. Three functions do the work, and none of them is in the repository:
- `InitializeHmmParameters()` sets the starting parameters;
- `EstimateHmmParameters()` estimates them (at most 50 iterations, tolerance $10^{-6}$);
- `ViterbiDecoding()` decodes the state path.

**Model quality.**

$$\mathrm{AIC} = -2\ell + 2(K^2 + 2K), \qquad \mathrm{BIC} = -2\ell + \ln(T)\,(K^2 + 2K),$$

Here $\ell$ is the log-likelihood and $T$ is the number of observations. Regime persistence is the mean of $\operatorname{diag}(P)$.

**Caveats:**
- The labels "Ultra-Accommodative", "Accommodative", "Neutral", "Restrictive" and "Crisis" are attached by state index. They are meaningful only if the estimation code orders the states accordingly, and this cannot be checked without that code.
- The "training history" animation is synthetic. Its frames interpolate linearly between a random row-stochastic matrix and the final estimate, and the plotted log-likelihood is the final value minus 2 per remaining frame. It does not show actual estimation iterations.
- This README reports no estimated transition probabilities or regime means, because the estimation code is not in the repository.

### Augmented Taylor rule, tails and Monte Carlo (`Main.R`; draft)

`Main.R` specifies the following steps. Only step 1 currently runs (see [Run order](#run-order)).

1. **Monthly panel.**
   - Monthly mean of RUONIA $r_t$, CPI inflation $\pi_t$, credit impulse, GDP level $Y_t$, the budget series, expected inflation and the MOEX + RTS volatility.
   - Missing values are filled with splines.
   - The GDP trend $Y_t^*$ is a LOESS fit (span 0.3), and the output gap is $g_t = 100\,(Y_t - Y_t^*)/Y_t^*$.
2. **Tails of the policy-rule residuals.** The script estimates the regression $r_t = b_0 + b_1 \pi_t + b_2 g_t + e_t$ by OLS. Then it fits:
   - an α-stable law to $\hat e_t$;
   - a Pareto law to the residuals above their 95th percentile;
   - the exceedance frequency $\hat\lambda$.
3. **Drift (augmented Taylor rule, after Taylor, 1993).** The planned regression is
   $$\Delta r_t = \beta_0 + \beta_1(\pi_t - \pi^*) + \beta_2 g_t + \beta_3 v_t + \beta_4 s_t + \beta_5 \hat E_t[\pi] + u_t.$$
   Here $\pi^*$ is the inflation target, $v_t$ is a volatility index, $s_t$ is a sanctions score and $\hat E_t[\pi]$ is estimated expected inflation. These inputs are not created yet.
4. **Shock parameters.** `stats4::mle()` maximises a sum of three log-likelihoods: a normal density, a symmetric α-stable density and a Poisson term for exceedances.
5. **Monte Carlo.** The script simulates 1,000 paths of 12 steps of
   $$r_t = r_{t-1} + \mu_t + \sigma\,\varepsilon_t + \sum_{j=1}^{N_t} J_{t,j}, \qquad \varepsilon_t \sim N(0,1), \quad N_t \sim \mathrm{Poisson}(\lambda), \quad J_{t,j} \sim S(\alpha, 0, \gamma, 0).$$
   Here $\mu_t$ is the fitted drift and $S(\alpha, 0, \gamma, 0)$ is a symmetric α-stable law with scale $\gamma$.
6. **Evaluation.** The script reports, at each horizon $h$, the empirical distribution $F_h$ of the simulated rates. The risk considered is a *rise* of the rate, so VaR and CVaR are taken in the upper tail. At confidence level $c$:

   $$\mathrm{VaR}_c(h) = F_h^{-1}(c) = \inf\{x : F_h(x) \ge c\}, \qquad \mathrm{CVaR}_c(h) = \mathbb{E}\left[r_{t+h} \mid r_{t+h} > \mathrm{VaR}_c(h)\right].$$

   - CVaR coincides with the conditional tail expectation (CTE) used in the article.
   - For downside risk, apply the same definitions to $-r_{t+h}$.
   - The script also computes the CRPS (Gneiting and Raftery, 2007), $\mathrm{CRPS}(F, y) = \int_{-\infty}^{\infty} \left(F(x) - \mathbf{1}\{x \ge y\}\right)^2 dx$, with `scoringRules::crps_sample()`.
   - It draws a fan chart of the 10/25/50/75/90 % quantiles.
   - As written, the code applies these statistics over the columns of the $12 \times 1000$ path matrix, that is, per path rather than per horizon. This is part of the unfinished draft.

### RDC feature selection (`2.Cpp/RDCScpp.cpp`)

The file exports two functions to R. It is inspired by the Randomized Dependence Coefficient (Lopez-Paz, Hennig and Schölkopf, 2013).

**`cancor_cpp(x, y)`** takes two numeric matrices with the same number of rows.
- It returns a zero result for degenerate input: fewer than 2 rows, no columns, or a column with variance below $10^{-12}$.
- Otherwise it returns `stats::cancor(x, y)` (`cor`, `xcoef`, `ycoef`, …).

**`RDCSelection_cpp(data, target = NULL, threshold = 0.1, k_features = 20, num_permutations = 500)`** takes `data`, a data frame of numeric columns. Each feature $x$ is compared with `target`; without a target, it is compared with a randomly chosen other feature. For each feature the function:
1. sets `NA` and infinite values to 0 and replaces values by ranks $R(\cdot)$;
2. builds $k$ = min(`k_features`, $n - 1$) random features $\sin(w\,R(x)/n)$ and $\cos(w\,R(y)/n)$, with a new $w \sim N(0, 1)$ for every observation and feature;
3. takes the largest canonical correlation $\hat\rho$;
4. computes the score:
   - with `num_permutations` $> 0$: $\hat\rho$ is recomputed on permutations of $y$, $f$ is the share of permutations whose correlation does not exceed $\hat\rho$, and the score is $|f - 0.5|/0.5$;
   - otherwise the score is $\hat\rho$.

Features whose score is at least `threshold` are kept. The function returns `data` (the selected columns plus the target), `rdc_scores`, `selected_features`, `removed_features` and `threshold`. The internal random generator has a fixed seed of 123, so the results are deterministic for given data.

How it differs from the original RDC:
- in the original, each random feature uses one projection weight (and bias) shared by all observations; here a new weight is drawn for every observation, and $x$ and $y$ use different non-linearities ($\sin$ and $\cos$);
- under independence, $f$ is roughly uniform, so the permutation score is roughly uniform on $[0, 1]$. A score close to 1 can also come from $f \approx 0$. Treat the score as a heuristic ranking, not as a test.

Example (run from the repository root):

```r
Rcpp::sourceCpp("2.Cpp/RDCScpp.cpp")

set.seed(1)
n  <- 300
df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = runif(n))
df$y <- sin(2 * df$x1) + 0.1 * rnorm(n)   # y depends on x1 only (non-linearly)

sel <- RDCSelection_cpp(df, target = "y", threshold = 0.1,
                        k_features = 20, num_permutations = 200)
sel$rdc_scores          # one score per feature, in [0, 1]
sel$selected_features   # features with score >= threshold, plus the target

cc <- cancor_cpp(as.matrix(df[, c("x1", "x2")]), as.matrix(df[, "y", drop = FALSE]))
cc$cor                  # canonical correlation(s), as returned by stats::cancor
```

## Roadmap (not in this repository yet)

The items below are planned or were described in earlier versions of this README. **None of them is implemented in this repository.**

- **Shared function modules** sourced by `8.SourceAll.R`: `HarVolatilityFunctions.R`, `HmmFunctions.R`, `BoostingFunctions.R` and `TalebR2Functions.R`. They include the HMM functions `InitializeHmmParameters()`, `EstimateHmmParameters()` and `ViterbiDecoding()`, and the SVD function `PredictTimeSeriesWithSVD()`.
- **Completing `Main.R`:**
  - an α-stable fitting step to replace `stable_fit()`;
  - the expectations module (`estimate_expectations()`);
  - the drift-regression inputs (`pi_target`, `volatility_index`, `sanctions_score`);
  - the Monte Carlo evaluation.
- **The factor layer from the article:** copula normalisation and SVD.
- **Gradient boosting** over candidate residual distributions (normal, Student-t, α-stable, GED) to choose the jump intensity, and a **Shapley decomposition** of R² for variable weights.
- **Rcpp acceleration** of Monte Carlo, HAR and HMM, and the helper functions announced earlier:
  - `OptimizedMonteCarlo()`, `FastHARCalculation()`, `FastHMMForwardBackward()`, `CompileRcppOptimizations()`;
  - `RunComprehensiveTests()`, `PrepareMainData()`, `AnalyzeFatTails()`, `IntegratedStochasticModel()`, `AddNewDataToModel()`;
  - `CreateMarkovNetworkVisualization()`, `CreateModelComparisonPlot()`, `CreateTransitionHeatmap()`.
- **More macro variables:** the exchange rate, commodity prices, fiscal and external balances, the labour market and financial conditions.
- **More methods:** VAR models, cointegration, Granger causality and structural-break tests.
- **Links to the Bank of Russia's quarterly projection models** of the Monetary Policy Department (ДДКП) and of the Research and Forecasting Department (ДИП), as a structural baseline and for calibration. The Bank of Russia publishes their code at https://www.cbr.ru/dkp/system_p/. It is not included here.

## References

Methods used in the code:

- Clements, A., & Preve, D. P. A. (2021). A practical guide to harnessing the HAR volatility model. *Journal of Banking & Finance*, 133, 106285. https://doi.org/10.1016/j.jbankfin.2021.106285 (cited in `Demo_Volatility_Analysis.R` as "Clements & Preve (2018)", the working-paper version).
- Corsi, F. (2009). A simple approximate long-memory model of realized volatility. *Journal of Financial Econometrics*, 7(2), 174–196. https://doi.org/10.1093/jjfinec/nbp001
- Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules, prediction, and estimation. *Journal of the American Statistical Association*, 102(477), 359–378. https://doi.org/10.1198/016214506000001437
- Hamilton, J. D. (1989). A new approach to the economic analysis of nonstationary time series and the business cycle. *Econometrica*, 57(2), 357–384. https://doi.org/10.2307/1912559
- Lopez-Paz, D., Hennig, P., & Schölkopf, B. (2013). The randomized dependence coefficient. *Advances in Neural Information Processing Systems 26*. https://proceedings.neurips.cc/paper/2013/hash/aab3238922bcc25a6f606eb525ffdc56-Abstract.html
- Taleb, N. N. (2020). *Statistical Consequences of Fat Tails: Real World Preasymptotics, Epistemology, and Applications*. STEM Academic Press. ISBN 978-1-5445-0805-4. https://arxiv.org/abs/2001.10488
- Taylor, J. B. (1993). Discretion versus policy rules in practice. *Carnegie-Rochester Conference Series on Public Policy*, 39, 195–214. https://doi.org/10.1016/0167-2231(93)90009-L

Background on tail risks in central banking:

- Cho, S., Mertens, T. M., & Williams, J. C. (2025, July 7). The zero lower bound remains a medium-term risk. *Liberty Street Economics*, Federal Reserve Bank of New York. https://libertystreeteconomics.newyorkfed.org/2025/07/the-zero-lower-bound-remains-a-medium-term-risk. It puts the seven-year-ahead ZLB risk at about 9 % on 27 May 2025; an earlier version of this README cited this figure as "NY Fed (July 2024)".
- Grishchenko, O. V., & Wilcox, L. (2024). Tale about inflation tails. *Finance and Economics Discussion Series* 2024-028, Board of Governors of the Federal Reserve System. https://doi.org/10.17016/FEDS.2024.028
- Kozlowski, J. (2024, July 2). Beyond the mean: Exploring tail risks in inflation expectations. *On the Economy*, Federal Reserve Bank of St. Louis. https://www.stlouisfed.org/on-the-economy/2024/jul/exploring-tail-risks-inflation-expectations
- BIS research on heterogeneous-agent New Keynesian models with ZLB constraints and on neural-network solution methods for DSGE models. The earlier README cited no specific publication.

## License

See [LICENSE](LICENSE).
