# CBModelsStohastic

[Русская версия](README.ru.md)

Research code in R and C++ for fat-tailed stochastic extensions of interest-rate modelling for Russia. It covers tail diagnostics, a fat-tail adjustment of R², HAR volatility, hidden Markov regimes and an augmented Taylor rule with Monte Carlo scenarios for the RUONIA rate.

Companion code for the article [*Stochastic modelling of monetary policy: a modular architecture for analysing extreme scenarios*](http://scipro.ru/article/07-03-2026) (Burakov, International Journal of Professional Science, 2026, No. 3-1, pp. 47–66).

## What is in this repository

The article describes a modular pipeline:
- latent factors from copula normalisation and SVD;
- HAR / stochastic volatility with calibrated tail distributions;
- a hidden Markov model (HMM) for regime switching;
- an augmented Taylor rule as the policy reaction operator;
- Monte Carlo fan charts with tail-risk metrics (VaR / CTE);
- evaluation with Taleb's R².

| Component | Code | Status |
|---|---|---|
| Integrated model on the real data: SVD forecasts of the macro panel, HAR volatility, 7-state HMM, calibrated and extended Taylor rules, Monte Carlo with α-stable errors, VaR / CTE per horizon | [`1.R/Main.R`](1.R/Main.R), functions in [`1.R/1.functions/all_fun.R`](1.R/1.functions/all_fun.R) | Runs end to end on [`3.Data/`](3.Data) (about 35 s) |
| Shared functions (data loading, SVD, HAR, HMM, calibration, Taleb R², Shapley R², tail metrics, plots) | [`1.R/1.functions/all_fun.R`](1.R/1.functions/all_fun.R), loaded by [`1.R/8.SourceAll.R`](1.R/8.SourceAll.R) | Library |
| Regressions of RUONIA and CPI inflation on the other panel variables | [`1.R/DataAnalysisAndRegression.R`](1.R/DataAnalysisAndRegression.R) | Runs on the real data |
| HMM regimes of the key rate (5 states) | [`1.R/5.HiddenMarkingsModels.R`](1.R/5.HiddenMarkingsModels.R) | Runs on the real data |
| HAR volatility of RUONIA with a stochastic band | [`1.R/Demo_HAR_Volatility.R`](1.R/Demo_HAR_Volatility.R) | Runs on the real data |
| SVD forecasting | [`1.R/Demo_SVD_Forecasting.R`](1.R/Demo_SVD_Forecasting.R) | Runs on simulated data |
| Tail diagnostics of distributions | [`1.R/3.dist_bust.R`](1.R/3.dist_bust.R), [`1.R/10.DashboardData.R`](1.R/10.DashboardData.R) | Run on simulated data |
| Fat-tail adjustment of R² (Taleb) | [`1.R/4.LmR2BoostingStohastic.R`](1.R/4.LmR2BoostingStohastic.R) | Runs on simulated data |
| HAR volatility, volatility of volatility, GARCH(1,1) benchmark | [`1.R/Demo_Volatility_Analysis.R`](1.R/Demo_Volatility_Analysis.R) | Runs on simulated data |
| RDC-based feature selection (C++) | [`2.Cpp/RDCScpp.cpp`](2.Cpp/RDCScpp.cpp) | Compiles with Rcpp; no R script calls it |
| Monte Carlo, HAR and HMM helpers (C++) | [`2.Cpp/FastMonteCarlo.cpp`](2.Cpp/FastMonteCarlo.cpp) | Compiled by `all_fun.R` when a C++ toolchain is available |
| Copula normalisation | — | Not in the repository: the SVD step standardises with z-scores |

All scripts listed above finish without errors when run with `Rscript` from the repository root (tested with R 4.4.1 on macOS).

## Corrections after publication

The repository version includes corrections made after the article was published. Numbers in the article that were computed with the earlier code can therefore differ from what this version prints.

- **Taleb's R².** `CalculateRobustR2()` (in `all_fun.R` and `4.LmR2BoostingStohastic.R`) applied the factor $\nu/(\nu-2)$ twice. It is now applied once (see [Taleb's R²](#fat-tail-adjustment-of-r-4lmr2boostingstohasticr)). `Main.R` now prints a Taleb R² of 0.89 for the base Taylor rule and 0.97 for the extended one, instead of 0.78 and 0.94. The Shapley weights and the Taleb R² values in `DataAnalysisAndRegression.R` change as well.
- **HMM log-likelihood.** `ForwardBackward()` returned only the log of the first normalising constant. This affected the EM stopping rule and the AIC / BIC in `5.HiddenMarkingsModels.R`. In `Main.R` the estimated HMM parameters are unchanged to three decimals. The simulated quantiles are unchanged to six decimals.
- **Tail risk and calibration.** VaR / CTE are computed per horizon across paths, and CRPS / PIT functions were added. The earlier draft computed these statistics per path and called `pit_sample()`, which does not exist.
- **RDC feature selection.** The score in `RDCSelection_cpp()` was replaced by the RDC and a permutation p-value.
- **α-stable sampler in `FastMonteCarlo.cpp`.** The Chambers–Mallows–Stuck formula was wrong, and every draw was positive. `Main.R` does not use this sampler.

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

- **R ≥ 4.0.** All scripts were tested with R 4.4.1 on macOS.
- **R packages.** [`1.R/9.libraries.r`](1.R/9.libraries.r) installs any missing package from CRAN and attaches all of them. `all_fun.R` (and therefore `8.SourceAll.R`) sources it.
  - `3.dist_bust.R` sources `9.libraries.r` itself. `4.LmR2BoostingStohastic.R`, `10.DashboardData.R` and `Demo_Volatility_Analysis.R` attach their own packages with `library()`.
  - `EstimateHarModel()` also uses `car` and `lmtest`, and the optional XGBoost comparison in `DataAnalysisAndRegression.R` runs only if `xgboost` is installed.
  - `qs` is optional. CRAN archived it on 2026-01-17. The scripts save their results with `SaveProjectResults()`: it calls `qs::qsave()` (`.qs` file) if `qs` is installed, and `saveRDS()` (`.rds` file with the same name) otherwise. `BuildDistDashboard(save_qs_path = …)` still needs `qs`.
- **C++ toolchain** (for [`2.Cpp/`](2.Cpp)). `Rcpp::sourceCpp()` needs a C++ compiler:
  - macOS: Xcode Command Line Tools (`xcode-select --install`);
  - Windows: the Rtools version that matches your R version;
  - Linux: `g++` and the R development headers (for example, `build-essential` and `r-base-dev` on Debian/Ubuntu).
  
  `all_fun.R` compiles `FastMonteCarlo.cpp` when it is sourced. If compilation fails, it prints a message and the R implementations are used.
- **Working directory.** All file paths in the scripts are relative to the repository root. Open [`CB_models_stohastic.Rproj`](CB_models_stohastic.Rproj) in RStudio, call `setwd()` with the repository root, or run `Rscript 1.R/<script>.R` from the root.
- **Locale.** Three data files have Cyrillic names, and the scripts read them by literal path. Run R in a UTF-8 locale (for example, `LANG=en_US.UTF-8`). In the C locale, or in R before 4.2 on Windows, these paths do not resolve.

## Data

[`3.Data/`](3.Data) contains the input files. `prep_data()` gives https://www.cbr.ru/dkp/statistics/ as the source of its inputs. The article (section "Materials and methods") lists these sources for its monthly panel:
- Bank of Russia: RUONIA and the key rate;
- Rosstat and Bank of Russia surveys: inflation and inflation expectations;
- Rosstat and the Ministry of Finance: GDP and the credit impulse;
- Moscow Exchange: market volatility.

The repository does not record the source of each file, except where a file names its publisher in its own header.

| File | Content | Frequency | Period | Used by |
|---|---|---|---|---|
| `RC_F11_01_2010_T04_07_2025.xlsx` | RUONIA: rate (`ruo`), volume (`vol`) and related fields | Daily | 11 Jan 2010 – 12 Aug 2025 | `prep_data()` (monthly mean of `ruo`), `Main.R` (realised values for CRPS / PIT), `Demo_HAR_Volatility.R` |
| `Инфляция и ключевая ставка Банка России_F01_12_2013_T18_06_2025.xlsx` | Key rate (% p.a.), CPI inflation (% y/y), inflation target | Monthly | Dec 2013 – Aug 2025 | `prep_data()`, `5.HiddenMarkingsModels.R`, `Demo_HAR_Volatility.R` |
| `add_info.xlsx` | Bank of Russia enterprise monitoring (30 sheets). Sheet 2: the Bank of Russia business-climate indicator (total, current, expected) and seven balance-of-responses series (output, demand, prices, costs) | Monthly | Sheet 2: Jan 2002 – Apr 2025 | `prep_data()` (sheet 2, renamed `B001`…`B0010`) |
| `bal_of_payments_ap.xlsx` | Balance of payments of the Russian Federation, analytical presentation: 52 indicators (`ind`) × quarters | Quarterly | Q1 1994 – Q1 2025 | `prep_data()` (renamed `BP001`…; each quarter is dated to its third month) |
| `кредитный_импульс_помесячно.csv` | Credit impulse | Monthly | Feb 2010 – Dec 2024 | `prep_data()`, `5.HiddenMarkingsModels.R` |
| `vvp_pomesyachno.csv` | GDP level (`ВВП_месяц`) | Monthly | Jan 2011 – Mar 2025 | `prep_data()`, `5.HiddenMarkingsModels.R` |
| `ожидаемая_инфляция.csv` | Expected inflation, % | Irregular until Mar 2020, monthly afterwards | Dec 2009 – Jun 2025 | `prep_data()`, `5.HiddenMarkingsModels.R` |
| `budget.xlsx` | Budget series (`Value`), called the "budget impulse" in the code | Monthly | Jan 2010 – Mar 2025 | `prep_data()` |
| `assets.xlsx` | MOEX and RTS index levels | Daily | 3 Jan 2017 – 30 Aug 2025 | `prep_data()` (volatility, see below) |
| `sunctions.csv` | Sanctions indicator: unnamed columns, non-decreasing values | Irregular | 21 Feb 2022 – 15 Jan 2025 | Read by `prep_data()`, but not merged (the merge line is commented out) |
| `stavka.xlsx` | RUONIA export with the same layout as `RC_F11_…` | Daily | 1 Nov 2010 – 16 Jun 2025 | Unused |
| `VVP_kvartal.xlsx` | GDP tables: current prices, constant prices, seasonally adjusted, volume indices, deflators | Quarterly | 1995 – 2025 | Unused |
| `bb.xlsx` | Bank of Russia Statistics Department: financial condition statistics of non-financial organisations | — | Header dated Dec 2024 | Unused |
| `exchange_rate.xlsx` | Main derived indicators of rouble exchange-rate dynamics | Monthly, quarterly | Jan 2005 – May 2025 | Unused |
| `inflation.xlsx` | CPI growth, % m/m, seasonally adjusted, by component | Monthly | From Jan 2002 | Unused |

**Monthly panel (`prep_data()`).**
- The panel covers Feb 2017 – Aug 2025 (103 months). It is anchored on the market volatility series: for each day, the mean absolute deviation of MOEX + RTS (the sum of the two index levels) over the last 30 trading days, times $\sqrt{30}$, averaged by month.
- The other files are merged by month. The August 2025 RUONIA value is the mean over 1–12 August, because the file ends on 12 August.
- Every missing value is replaced by the fitted value of a linear regression of that series on time over the whole panel. Series that end earlier, such as GDP after March 2025, are therefore extended by a linear trend.

Scripts write their results to `3.Data/*.qs` or `3.Data/*.rds`. These files are listed in `.gitignore`.

## Run order

From the repository root:

```bash
Rscript 1.R/Main.R                       # integrated model; saves 3.Data/final_integrated_cb_model.qs (.rds)
Rscript 1.R/5.HiddenMarkingsModels.R     # 5-state HMM of the key rate
Rscript 1.R/DataAnalysisAndRegression.R  # regressions of RUONIA and CPI on the panel
Rscript 1.R/Demo_HAR_Volatility.R        # HAR volatility of RUONIA with a stochastic band
Rscript 1.R/Demo_SVD_Forecasting.R       # SVD forecasting on simulated data
Rscript 1.R/3.dist_bust.R                # tail diagnostics on simulated data
Rscript 1.R/4.LmR2BoostingStohastic.R    # standard vs fat-tail-adjusted R² on simulated data
Rscript 1.R/10.DashboardData.R           # distribution dashboard on six simulated samples
Rscript 1.R/Demo_Volatility_Analysis.R   # HAR vs GARCH(1,1) on simulated returns
```

Plots are created as plotly objects. To see them, run the scripts in an interactive session, for example `source("1.R/Main.R", print.eval = TRUE)`, and then print `beautiful_plot`, `model_comparison_plot`, `transition_heatmap`, and so on.

The C++ files can be compiled on their own:

```r
Rcpp::sourceCpp("2.Cpp/RDCScpp.cpp")        # cancor_cpp(), RDCSelection_cpp()
Rcpp::sourceCpp("2.Cpp/FastMonteCarlo.cpp") # Fast*() helpers used by all_fun.R
```

## Scripts

```
CBModelsStohastic/
├── CB_models_stohastic.Rproj        RStudio project; opening it sets the working directory to the repository root
├── README.md, README.ru.md          Documentation (English, Russian)
├── 1.R/
│   ├── 1.functions/all_fun.R        Shared functions (see below); sources 9.libraries.r and compiles 2.Cpp/FastMonteCarlo.cpp
│   ├── 3.dist_bust.R                Tail diagnostics: normal vs stochastic-volatility vs mixture samples
│   ├── 4.LmR2BoostingStohastic.R    Standard vs fat-tail-adjusted R² (Taleb) on simulated regressions
│   ├── 5.HiddenMarkingsModels.R     5-regime HMM for the key rate
│   ├── 8.SourceAll.R                Sources all_fun.R and reports a missing file or an error
│   ├── 9.libraries.r                Installs and attaches all R packages
│   ├── 10.DashboardData.R           BuildDistDashboard(): moments, normality tests, Hill / power-law tails
│   ├── DataAnalysisAndRegression.R  OLS of RUONIA and CPI on the panel, Taleb R², optional XGBoost comparison
│   ├── Demo_HAR_Volatility.R        HAR forecast of RUONIA volatility as in Main.R, with a stochastic band
│   ├── Demo_SVD_Forecasting.R       PredictTimeSeriesWithSVD() on simulated macro series
│   ├── Demo_Volatility_Analysis.R   HAR volatility with shocks and vol-of-vol band; GARCH(1,1) benchmark
│   └── Main.R                       Integrated model (see "Integrated model" below)
├── 2.Cpp/
│   ├── FastMonteCarlo.cpp           α-stable and compound-Poisson draws, integrated simulation, HAR windows, HMM forward-backward
│   └── RDCScpp.cpp                  cancor_cpp() and RDCSelection_cpp(): RDC-based feature selection
└── 3.Data/                          Input data (see Data)
```

The file name `5.HiddenMarkingsModels.R` is kept as is: it means "Hidden Markov Models".

Functions in `1.R/1.functions/all_fun.R`:

| Group | Functions |
|---|---|
| Data | `prep_data()` |
| SVD forecasting | `PredictTimeSeriesWithSVD()`, `ExtrapolateLatentFactorWithSeasonality()`, `InterpolateDataWithSVD()`, `SmartExtrapolate()` |
| HAR | `EstimateHarModel()`, `ForecastHarVolatility()`, `EstimateStableHarModel()`, `ForecastStableHarVolatility()` |
| HMM | `InitializeHmmParameters()`, `ForwardAlgorithm()`, `BackwardAlgorithm()`, `ForwardBackward()`, `EstimateHmmParameters()`, `ViterbiDecoding()`, `PredictRegimeProbabilities()` |
| Calibration | `CalibrateStochasticDistributions()`, `CalibrateTaylorRuleParameters()`, `EstimateAlphaParameter()`, `SelectBestDistribution()` |
| R² | `CalculateRobustR2()`, `EstimateStudentScale()`, `ComputeShapleyR2()` |
| Tails and evaluation | `AnalyzeFatTails()`, `ComputeTailRiskByHorizon()`, `EvaluatePredictiveSample()` |
| Simulation | `IntegratedStochasticCBModel()` (`Main.R` defines its own version), `OptimizedMonteCarlo()`, `FastHARCalculation()`, `FastHMMForwardBackward()`, `CompileRcppOptimizations()` |
| Plots | `CreateBeautifulMonteCarloPlot()`, `CreateMacroTrendsPlot()` |
| Tests and utilities | `TestAllFunctions()`, `RunComprehensiveTests()`, `SaveProjectResults()`, `SaveResults()`, `LoadResults()` |

`Main.R` defines the plotting and comparison functions `CreateBaselineForecasts()`, `CreateModelComparisonPlot()`, `CreateMarkovNetworkVisualization()`, `CreateTransitionHeatmap()`, `CreateRegimeEvolutionPlot()`, `CreateMacroForecastPlots()` and `CreateAdvancedMarkovChainPlot()`.

## Method

This section describes what the code implements. Symbols are used consistently throughout:

| Symbol | Meaning |
|---|---|
| $\nu$ | Tail index / Student-t degrees of freedom (`ALPHA_TAIL`, `alpha_tail` in the code) |
| $\alpha$ | Stability index of an α-stable law |
| $c$ | Confidence level of VaR / CTE (`Main.R` uses $c = 0.95$) |
| $K$ | Number of HMM states |
| $P$ | HMM transition matrix |
| $h$ | Forecast horizon in months |

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

**Taleb's formula.** For $y = a x + s\,\varepsilon$ with $\operatorname{Var}(x) = 1$ and a standard Student-t error $\varepsilon \sim t_\nu$, the error variance is $s^2\nu/(\nu-2)$. The population R² is then

$$R^2 = \frac{a^2}{a^2 + s^2\,\nu/(\nu-2)} = \frac{\tilde a^2(\nu-2)}{\tilde a^2(\nu-2) + \nu}, \qquad \tilde a^2 = a^2/s^2,$$

which is the expression Taleb (2020) gives for $s = 1$.

**R² variants.** After OLS, `CalculateRobustR2(y, ŷ, ν)` computes three values. Here $\hat\varepsilon_i = y_i - \hat y_i$, and $\hat s$ is the scale of a Student-t($\nu$) fit to $\hat\varepsilon$ with $\nu$ fixed (`EstimateStudentScale()`: maximum likelihood with `MASS::fitdistr()`; $\operatorname{MAD}(\hat\varepsilon)/t_\nu^{-1}(0.75)$ if that fit fails):

$$R^2_{\text{std}} = 1 - \frac{\sum_i \hat\varepsilon_i^2}{\sum_i (y_i - \bar y)^2},$$

$$R^2_{\text{Taleb}} = \frac{\widehat{\operatorname{Var}}(\hat y)}{\widehat{\operatorname{Var}}(\hat y) + \frac{\nu}{\nu - 2}\,\hat s^2}, \qquad R^2_{\text{MLE}} = \max\left(0,\; 1 - \frac{n\,\frac{\nu}{\nu - 2}\,\hat s^2}{\sum_i (y_i - \bar y)^2}\right).$$

- The factor $\nu/(\nu-2)$ converts the squared t scale into a variance, so it is applied to $\hat s^2$ once. The sample variance of $\hat\varepsilon$ already estimates $s^2\nu/(\nu-2)$ and is not multiplied by the factor again. Earlier versions of the code did multiply it again (see [Corrections after publication](#corrections-after-publication)).
- $R^2_{\text{MLE}}$ replaces the residual mean square $\sum_i\hat\varepsilon_i^2/n$ in $R^2_{\text{std}}$ by the error variance of the fitted t model.
- Both adjusted values are set to 0 when $\nu \le 2.1$, and $R^2_{\text{Taleb}}$ is clipped to $[0, 1]$. The `all_fun.R` version instead uses $\max(\nu, 2.1)$ in $R^2_{\text{Taleb}}$, returns $R^2_{\text{std}}$ as $R^2_{\text{MLE}}$ when $\nu \le 2$, and clips all three values to $[0, 1]$.
- The script uses $\nu = 30$ for normal residuals, $\nu = 2.8$ for fat-tailed residuals and $\nu = \max(4,\, 0.8 \cdot 2.8) = 4$ for the mixture.

**Behaviour.**
- For Gaussian errors with $\nu = 1000$, $R^2_{\text{Taleb}}$ equals the classic R² (0.4972 vs 0.4972 in a simulation with $n = 10\,000$).
- For Student-t errors with $\nu = 3, 4, 6$, it is within 0.01 of the population R².
- The script's fat-tailed and mixture residuals are not Student-t: the 4 % outliers and the mixture inflate the sample variance, but the t fit is driven by the bulk of the data. For these two residual types the script prints adjusted values above the standard R² (0.82 vs 0.70 and 0.90 vs 0.81).

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

**Stochastic band.** The script simulates 1,000 paths of a mean-reverting volatility process around the HAR forecast. The code calls it Heston-type and attributes it to Lu (2022), for which it gives no full reference. The parameters are $\kappa = 0.5$, $\sigma_v = 0.3$, leverage $-0.5$, and the weight of the HAR forecast is $0.6 + 0.2e^{-t/5}$. It reports quantiles from 1 % to 99 % for each day, across paths.

**Benchmark and accuracy.** The benchmark is a GARCH(1,1) model with normal innovations (`rugarch`), annualised as $\hat\sigma_t\sqrt{252}$. The script compares MSE, MAE and MAPE on the test window.

### SVD forecasting (`PredictTimeSeriesWithSVD()`; `Demo_SVD_Forecasting.R`, `Main.R`)

For a table of numeric series:
1. Each series with at least two seasonal periods is split by STL (`stl(s.window = "periodic", robust = TRUE)`) into trend and seasonal parts.
2. Sine and cosine terms with the annual, semi-annual and quarterly frequency are added. The matrix is standardised with z-scores (no copula transform).
3. The SVD keeps the smallest rank whose singular values cover at least 95 % of the energy $\sum_j d_j^2$.
4. Each retained left singular vector is extrapolated by `ExtrapolateLatentFactorWithSeasonality()`.
   - The function first tries `vars::VAR()` on the component. `vars::VAR()` rejects a single series ("The matrix 'y' should contain at least two variables"), so in practice the error handler takes over.
   - The handler fits `ar(method = "yule-walker", order.max = 2)`. If that also fails, it falls back to `forecast::ets()` and then to a linear trend.
5. The components are recombined through the right singular vectors and de-standardised. Each series is returned as trend + seasonal part.

`Demo_SVD_Forecasting.R` applies it to five simulated monthly series (120 months; the last 24 are held out). `Main.R` applies it to the 70 numeric columns of the panel with a horizon of 19 months (rank 8, 95.8 % of the energy on the current data).

### HMM regimes (`5.HiddenMarkingsModels.R`)

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

**Model.** A Markov-switching model in the sense of Hamilton (1989) with $K = 5$ hidden states $S_t$, transition matrix $P_{jk} = P(S_{t+1} = k \mid S_t = j)$ and Gaussian emissions $r_t \mid S_t = k \sim N(\mu_k, \sigma_k^2)$. The functions are in `all_fun.R`:
- `InitializeHmmParameters()`: k-means starting values, $P$ with 0.8 on the diagonal;
- `EstimateHmmParameters()`: Baum–Welch EM (here at most 50 iterations, tolerance $10^{-6}$), with $\sigma_k \ge 0.1$. After estimation the states are sorted by $\mu_k$ in ascending order;
- `ViterbiDecoding()`: the most likely state path.

**Model quality.**

$$\mathrm{AIC} = -2\ell + 2(K^2 + 2K), \qquad \mathrm{BIC} = -2\ell + \ln(T)\,(K^2 + 2K),$$

Here $\ell$ is the log-likelihood, the sum of the log normalising constants of the scaled forward pass, and $T$ is the number of observations. Regime persistence is the mean of $\operatorname{diag}(P)$.

**Caveats:**
- The labels "Ultra-Accommodative", "Accommodative", "Neutral", "Restrictive" and "Crisis" are attached by state index. Because the states are sorted by $\mu_k$, the labels follow the mean key rate of each state from lowest to highest.
- The "training history" animation is synthetic. Its frames interpolate linearly between a random row-stochastic matrix and the final estimate, and the plotted log-likelihood is the final value minus 2 per remaining frame. It does not show actual estimation iterations.

### Integrated model (`Main.R`)

All constants below are set in `Main.R` or `all_fun.R`. Horizon $H = 19$ months, $N = 10\,000$ paths, inflation target $\pi^* = 4\,\%$.

1. **Panel and forecasts of the drivers.** `prep_data()` builds the panel (see [Data](#data)). The training sample ends at the last month with RUONIA, CPI and GDP (Aug 2025 on the current data, after the gap filling). `PredictTimeSeriesWithSVD()` forecasts all numeric columns for $H$ months.
   - Derived variables: $\pi_t - \pi^*$ (`inflation_deviation`); the output gap $g_t = 100\,(Y_t - Y_t^*)/Y_t^*$ with a LOESS trend $Y_t^*$ (span 0.3); the 6-month rolling standard deviation of RUONIA $v_t$ (`rate_volatility`); credit impulse / GDP; the key rate minus inflation.
2. **Distribution calibration** (`CalibrateStochasticDistributions()`), by maximum likelihood:
   - volatility: $\Delta v_t = \kappa(\theta - v_{t-1}) + \sigma_v\sqrt{v_{t-1}}\,\eta_t$, $\eta_t \sim N(0,1)$;
   - CPI inflation and the output gap: Student-t or Pareto, chosen by AIC;
   - residuals of the OLS regression of RUONIA on $\pi_t - \pi^*$ and $g_t$: an α-stable law $S(\alpha, \beta, \gamma, \delta)$ (`stabledist::dstable`).
3. **HAR volatility.** `highfrequency::HARmodel(periods = c(5, 22))` on the monthly $v_t$. The model is re-estimated after each step, and 19 one-step forecasts are chained.
4. **HMM regimes.**
   - Observations: a composite indicator, the weighted sum of the z-scores of all numeric panel columns. `gdp_gap`, `B002`, `credit_to_gdp`, `budget_impulse`, `policy_tightness` and `ruonia_rate` have weight 0.07 each, and the remaining columns share 0.58 equally.
   - $K = 7$ states, Baum–Welch EM with at most 1,000 iterations, then Viterbi decoding.
   - Forecast regime probabilities start from state 7: $p(h) = e_7^\top P^{\,h-1}$ (`PredictRegimeProbabilities(current_state = 7)`).
5. **Taylor rule.**
   - `CalibrateTaylorRuleParameters()` fits $i_t = \rho\, i_{t-1} + (1-\rho)\,[r^* + \phi_\pi(\pi_t - \pi^*) + \phi_y g_t] + \varepsilon_t$, $\varepsilon_t \sim N(0, \sigma^2)$, by maximum likelihood.
   - `Main.R` then keeps $\rho = 0.8$ instead of the fitted value and lowers $r^*$ by 2 percentage points (`NEUTRAL_RATE_ADJUSTMENT`).
   - The extended rule is an OLS regression of RUONIA on its lag, $\pi_t - \pi^*$, $g_t$, the HAR volatility, the regime-7 dummy, expected inflation and credit / GDP.
   - The Hill estimate of the tail index of the base-rule residuals, clipped to $[1.1, 4]$ (`EstimateAlphaParameter()`), is used as $\nu$ in Taleb's R² and as $\alpha$ of the simulated errors.
   - `ComputeShapleyR2()` gives a simplified Shapley split of Taleb's R². For each variable it averages its R² alone with the drop in R² when it is removed from the full model, and then normalises the values to sum to 1.
6. **Monte Carlo** (`IntegratedStochasticCBModel()` as defined in `Main.R`). For $t = 1, \dots, H$ on each path:
   - volatility: $V_t = 0.6\,\mathrm{HAR}_t + 0.4\,[V_{t-1} + \kappa(\theta - V_{t-1}) + \sigma_v\sqrt{V_{t-1}}\,\eta_t]$;
   - drivers: $\tilde\pi_t = 0.7\,\hat\pi_t + 0.3\,\xi^\pi_t$ and $\tilde g_t = 0.7\,\hat g_t + 0.3\,\xi^g_t$, where $\hat\pi_t, \hat g_t$ are the SVD forecasts and $\xi$ are draws from the calibrated t or Pareto laws;
   - neutral rate: $r^*_t = r^*_0 + (7 - r^*_0)(1 - e^{-0.35t}) + \sum_k \mu_k\, p_k(t)$, where $\mu_k$ are the HMM emission means and $p_k(t)$ the forecast regime probabilities;
   - $L_t = r^*_t + \sum_j b_j e^{-\kappa_j t} X_{j,t}$: the extended-rule coefficients $b_j$ of $\tilde\pi_t$, $\tilde g_t$, $V_t$, expected inflation $-\,\pi^*$ and credit / GDP, damped with $\kappa_j = 0.25, 0.15, 0.20, 0.20, 0.10$;
   - $r_t = \min\big(40, \max(0.1,\; \rho_{\text{ext}}\, r_{t-1} + (1-\rho_{\text{ext}})\,L_t + e_t)\big)$, with $\rho_{\text{ext}}$ = the lag coefficient of the extended rule clipped to $[0.4, 0.85]$, and $e_t \sim S(\alpha, 0, 0.5\sqrt{V_t}, 0)$ with $\alpha$ clipped to $[1.01, 1.99]$.
   
   The simulated regime path starts from the last Viterbi state and follows $P$. It is stored, but the rate equation uses only $p_k(t)$. The α-stable parameters calibrated in step 2 are not used in the simulation.
7. **Evaluation.** For each horizon $h$, $F_h$ is the empirical distribution of the 10,000 simulated rates $r_{t+h}$. All statistics below are computed across paths for each horizon.
   - **Fan chart:** the 5, 10, 25, 50, 75, 90 and 95 % quantiles of $F_h$.
   - **VaR / CTE** (`ComputeTailRiskByHorizon()`). The risk considered is a *rise* of the rate, so both are taken in the upper tail. At confidence level $c$:

     $$\mathrm{VaR}_c(h) = F_h^{-1}(c) = \inf\{x : F_h(x) \ge c\}, \qquad \mathrm{CTE}_c(h) = \mathbb{E}\left[r_{t+h} \mid r_{t+h} > \mathrm{VaR}_c(h)\right].$$

     CTE is the same as CVaR here. For downside risk, apply the same definitions to $-r_{t+h}$.
   - **CRPS and PIT** (`EvaluatePredictiveSample()`), for horizons with a realised value $y_h$:
     - $\mathrm{CRPS}(F_h, y_h) = \int (F_h(x) - \mathbf 1\{x \ge y_h\})^2\,dx$ (Gneiting and Raftery, 2007), computed with `scoringRules::crps_sample()`;
     - $\mathrm{PIT}_h = F_h(y_h)$, randomised as $F_h(y_h^-) + U\,[F_h(y_h) - F_h(y_h^-)]$ only when simulated values tie with $y_h$ (Czado, Gneiting and Held, 2009).

     `scoringRules` has no PIT function, so PIT is computed directly. `Main.R` takes the realised values from the complete months of the RUONIA file. The current file ends before the forecast period, so this step is skipped with a message.
   - A deterministic benchmark (`CreateBaselineForecasts()`): 0.3 × the mean of the last 21 months + 0.4 × the linear trend of the last 24 months + 0.3 × a Taylor rule with decaying gaps.
   - The results are saved in `final_results` (`3.Data/final_integrated_cb_model.qs` or `.rds`), including `tail_risk` and `forecast_evaluation`.

### HAR volatility of RUONIA (`Demo_HAR_Volatility.R`)

The script repeats the volatility step of `Main.R` on a smaller panel:
- monthly RUONIA and the key-rate / inflation file from December 2013;
- $v_t$ is the 6-month rolling standard deviation of RUONIA;
- the HAR model is the same as in step 3 above, with 19 chained one-step forecasts;
- 1,000 paths of $V_t = 0.6\,\mathrm{HAR}_t + 0.4\,[V_{t-1} + 0.5(\theta - V_{t-1}) + 0.3\sqrt{V_{t-1}}\,\eta_t]$ give quantiles per horizon, with $\theta$ the sample mean of $v_t$.

### Regressions on the panel (`DataAnalysisAndRegression.R`)

For RUONIA and CPI inflation, the script runs an OLS regression on all "raw" numeric panel columns, that is, all except the derived ones listed in the script. It reports R², adjusted R², AIC, BIC and Taleb's R², with $\nu$ = the Hill estimate of the residual tail.
- On the current panel the number of regressors is close to the number of months, so `lm()` drops aliased regressors and `predict()` warns about a rank-deficient fit.
- A random 80/20 split compares the OLS test RMSE with XGBoost if `xgboost` is installed.

### RDC feature selection (`2.Cpp/RDCScpp.cpp`)

The file exports two functions to R. They follow the Randomized Dependence Coefficient (RDC) of Lopez-Paz, Hennig and Schölkopf (2013).

**`cancor_cpp(x, y)`** takes two numeric matrices with the same number of rows.
- It returns a zero result for degenerate input: fewer than 2 rows, no columns, or a column with variance below $10^{-12}$.
- Otherwise it returns `stats::cancor(x, y)` (`cor`, `xcoef`, `ycoef`, …).

**`RDCSelection_cpp(data, target = NULL, threshold = 0.1, k_features = 20, num_permutations = 500)`** takes `data`, a data frame of numeric columns. Each feature $x$ is compared with `target`; without a target, it is compared with a randomly chosen other feature. For each pair $(x, y)$ the function:
1. sets `NA` and infinite values to 0 and computes the copula transform $u = R(x)/n$ (ranks, ties averaged);
2. builds $k$ = min(`k_features`, $n - 1$) random features $\sin\big(\tfrac{s}{2}(w_i u + b_i)\big)$ with $(w_i, b_i) \sim N(0, 1)$ and $s = 1/6$. Each feature has one pair $(w_i, b_i)$, shared by all observations. The same is done for $y$ with its own weights;
3. takes the RDC, the largest canonical correlation $\hat\rho$ between the two feature sets (`rdc_scores`);
4. with `num_permutations` $= B > 0$, permutes the rows of the $y$ features $B$ times and returns the p-value $p = \big(1 + \#\{\hat\rho_{\text{perm}} \ge \hat\rho\}\big)/(1 + B)$ (`p_values`).

**Selection.** With permutations, a feature is kept when $p \le$ `threshold`, so `threshold` is a significance level. Independent features are then kept with probability about `threshold`. Without permutations, a feature is kept when $\hat\rho \ge$ `threshold`. The RDC of a finite sample is positive even under independence (about 0.2 for the independent features in the example below), so the permutation mode is the one to use for selection.

The function returns `data` (the selected columns plus the target), `rdc_scores`, `p_values` (`NA` without permutations), `selected_features`, `removed_features` and `threshold`. The internal random generator has a fixed seed of 123, so the results are deterministic for given data.

In a simulation with 50 data sets ($n = 200$, threshold 0.1, 100 permutations), a linear and a non-linear regressor were kept in all 50 data sets, and two independent regressors in 10 % and 12 % of them.

Example (run from the repository root):

```r
Rcpp::sourceCpp("2.Cpp/RDCScpp.cpp")

set.seed(1)
n  <- 300
df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = runif(n))
df$y <- sin(2 * df$x1) + 0.1 * rnorm(n)   # y depends on x1 only (non-linearly)

sel <- RDCSelection_cpp(df, target = "y", threshold = 0.1,
                        k_features = 20, num_permutations = 200)
sel$rdc_scores          # x1 0.975, x2 0.208, x3 0.193
sel$p_values            # x1 0.005, x2 0.229, x3 0.522
sel$selected_features   # "x1" "y"

cc <- cancor_cpp(as.matrix(df[, c("x1", "x2")]), as.matrix(df[, "y", drop = FALSE]))
cc$cor                  # canonical correlation(s), as returned by stats::cancor
```

### C++ helpers (`2.Cpp/FastMonteCarlo.cpp`)

- `FastAlphaStableSimulation()` draws $S(\alpha, \beta, \gamma, \delta)$ variables (parameterisation 1, as `stabledist::rstable(pm = 1)`). It uses the Chambers–Mallows–Stuck method in the form of Weron (1996).
- `FastCompoundPoissonSimulation()` draws compound Poisson sums of normal jumps.
- `FastIntegratedStochasticSimulation()` simulates paths with Markov regimes, weighted drivers, Gaussian shocks and Poisson-driven symmetric α-stable jumps. `OptimizedMonteCarlo()` calls it and returns a horizon × path matrix.
- `FastHarVolatilityCalculation()` computes 1/5/22/66-period moving averages.
- `FastForwardBackwardAlgorithm()` returns smoothed state probabilities, and `FastLogLikelihood()` the HMM log-likelihood.

These functions use `std::random_device` seeds, so their draws do not follow `set.seed()`. `Main.R` does not call them.

### Known limitations

- `TestAllFunctions()` and `RunComprehensiveTests()` call `PrepareHarData()`, which is not defined, so their HAR test always fails. Their error handlers assign to a local copy of the result list, so failed tests are missing from the returned list instead of being reported. The pass count can exceed the number of tests ("8 of 7").
- In `Main.R`, `CreateRegimeEvolutionPlot()` runs inside `tryCatch()`. On error it silently falls back to a simpler plot.
- On the current data, the α-stable fit in `CalibrateStochasticDistributions()` reaches the upper bound $\alpha = 2$.

## Roadmap (not in this repository yet)

The items below are planned or were described in earlier versions of this README. **None of them is implemented in this repository.**

- **Copula normalisation** of the factor layer described in the article. `PredictTimeSeriesWithSVD()` standardises with z-scores.
- **A multivariate VAR for the SVD factors.** The code tries `vars::VAR()` on each factor separately and falls back to a univariate AR model (see [SVD forecasting](#svd-forecasting-predicttimeserieswithsvd-demo_svd_forecastingr-mainr)).
- **The expectations module** `estimate_expectations()` (an RNN / Transformer estimate of $E_t[\pi]$). `Main.R` uses the survey series of expected inflation.
- **The jump-diffusion variant of the earlier draft of `Main.R`:** Poisson-driven α-stable jumps with $\sigma$, $\gamma$ and $\lambda$ fitted jointly by maximum likelihood, and a Pareto fit of the residual tail. The integrated `Main.R` uses α-stable errors without a jump component.
- **Gradient boosting** over candidate residual distributions. `SelectBestDistribution()` weights the normal, Student-t, α-stable and "GED" candidates by their likelihoods. The "GED" candidate is a normal law with a $\sqrt 2$ larger standard deviation, and no boosting is done.
- **The sanctions indicator** in the panel: `sunctions.csv` is read but not merged.
- **Helper functions announced earlier:** `PrepareHarData()`, `PrepareMainData()`, `IntegratedStochasticModel()` and `AddNewDataToModel()`.
- **More macro variables:** the exchange rate, commodity prices, fiscal and external balances, the labour market and financial conditions.
- **More methods:** cointegration, Granger causality and structural-break tests.
- **Links to the Bank of Russia's quarterly projection models** (КПМ) of the Monetary Policy Department (ДДКП) and of the Research and Forecasting Department (ДИП), as a structural baseline and for calibration. The Bank of Russia publishes their code (`kpm_ddkp.rar`, `kpm_dip.rar`) and information notes that describe the models at https://www.cbr.ru/dkp/system_p/. They are not included here.

## References

Methods used in the code:

- Chambers, J. M., Mallows, C. L., & Stuck, B. W. (1976). A method for simulating stable random variables. *Journal of the American Statistical Association*, 71(354), 340–344. https://doi.org/10.1080/01621459.1976.10480344
- Clements, A., & Preve, D. P. A. (2021). A practical guide to harnessing the HAR volatility model. *Journal of Banking & Finance*, 133, 106285. https://doi.org/10.1016/j.jbankfin.2021.106285 (cited in `Demo_Volatility_Analysis.R` as "Clements & Preve (2018)", the working-paper version).
- Corsi, F. (2009). A simple approximate long-memory model of realized volatility. *Journal of Financial Econometrics*, 7(2), 174–196. https://doi.org/10.1093/jjfinec/nbp001
- Czado, C., Gneiting, T., & Held, L. (2009). Predictive model assessment for count data. *Biometrics*, 65(4), 1254–1261. https://doi.org/10.1111/j.1541-0420.2009.01191.x
- Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules, prediction, and estimation. *Journal of the American Statistical Association*, 102(477), 359–378. https://doi.org/10.1198/016214506000001437
- Hamilton, J. D. (1989). A new approach to the economic analysis of nonstationary time series and the business cycle. *Econometrica*, 57(2), 357–384. https://doi.org/10.2307/1912559
- Lopez-Paz, D., Hennig, P., & Schölkopf, B. (2013). The randomized dependence coefficient. *Advances in Neural Information Processing Systems 26*. https://proceedings.neurips.cc/paper/2013/hash/aab3238922bcc25a6f606eb525ffdc56-Abstract.html
- Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. *Proceedings of the IEEE*, 77(2), 257–286. https://doi.org/10.1109/5.18626
- Taleb, N. N. (2020). *Statistical Consequences of Fat Tails: Real World Preasymptotics, Epistemology, and Applications*. STEM Academic Press. ISBN 978-1-5445-0805-4. https://arxiv.org/abs/2001.10488
- Taylor, J. B. (1993). Discretion versus policy rules in practice. *Carnegie-Rochester Conference Series on Public Policy*, 39, 195–214. https://doi.org/10.1016/0167-2231(93)90009-L
- Weron, R. (1996). On the Chambers–Mallows–Stuck method for simulating skewed stable random variables. *Statistics & Probability Letters*, 28(2), 165–171. https://doi.org/10.1016/0167-7152(95)00113-1

Background on tail risks in central banking:

- Cho, S., Mertens, T. M., & Williams, J. C. (2025, July 7). The zero lower bound remains a medium-term risk. *Liberty Street Economics*, Federal Reserve Bank of New York. https://libertystreeteconomics.newyorkfed.org/2025/07/the-zero-lower-bound-remains-a-medium-term-risk. It puts the seven-year-ahead ZLB risk at about 9 % on 27 May 2025; an earlier version of this README cited this figure as "NY Fed (July 2024)".
- Grishchenko, O. V., & Wilcox, L. (2024). Tale about inflation tails. *Finance and Economics Discussion Series* 2024-028, Board of Governors of the Federal Reserve System. https://doi.org/10.17016/FEDS.2024.028
- Kozlowski, J. (2024, July 2). Beyond the mean: Exploring tail risks in inflation expectations. *On the Economy*, Federal Reserve Bank of St. Louis. https://www.stlouisfed.org/on-the-economy/2024/jul/exploring-tail-risks-inflation-expectations
- BIS research on heterogeneous-agent New Keynesian models with ZLB constraints and on neural-network solution methods for DSGE models. The earlier README cited no specific publication.

## License

See [LICENSE](LICENSE).
