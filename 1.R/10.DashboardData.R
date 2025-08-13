# Distribution Dashboard: plots, moments, normality tests and power law tails
BuildDistDashboard <- function(x, nbins = 100, sample_shapiro = 5000,
                               tail_fraction = 0.1, conf_level = 0.95, save_qs_path = NULL) {
  suppressWarnings({
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Package data.table required")
    if (!requireNamespace("plotly", quietly = TRUE)) stop("Package plotly required")
    if (!requireNamespace("stringr", quietly = TRUE)) stop("Package stringr required")
  })
  # Need attach for := operations and interactive plots
  if (!is.null(save_qs_path)) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Package qs required for QS saving")
  }

  # Data preparation
  x_num <- suppressWarnings(as.numeric(x))
  x_num <- x_num[is.finite(x_num)]
  n <- length(x_num)
  if (n < 10) stop("Too few observations (<10)")

  # Basic statistics (without external packages)
  mean_x <- mean(x_num)
  sd_x <- sd(x_num)
  var_x <- var(x_num)
  med_x <- median(x_num)
  mad_x <- stats::mad(x_num, constant = 1.4826)
  qtls <- stats::quantile(x_num, probs = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99), names = FALSE, type = 7)
  z <- (x_num - mean_x) / sd_x
  skew_ex <- mean(z^3)
  kurt_ex <- mean(z^4) - 3
  jb_stat <- n/6 * (skew_ex^2 + (kurt_ex^2)/4)
  jb_p <- stats::pchisq(jb_stat, df = 2, lower.tail = FALSE)

  # Normality tests
  w_sample <- if (n > sample_shapiro) sample(x_num, sample_shapiro) else x_num
  shap <- tryCatch(stats::shapiro.test(w_sample), error = function(e) NULL)
  shap_p <- if (!is.null(shap)) shap$p.value else NA_real_
  ad_p <- tryCatch({
    if (requireNamespace("nortest", quietly = TRUE)) nortest::ad.test(x_num)$p.value else NA_real_
  }, error = function(e) NA_real_)

  # Tail characteristics (right tail)
  # Threshold for tail
  tail_thr <- stats::quantile(x_num, probs = 1 - tail_fraction, names = FALSE)
  x_tail <- x_num[x_num >= tail_thr]
  k_tail <- length(x_tail)
  # Hill estimator (alpha = 1 / gamma)
  hill_alpha <- NA_real_
  hill_curve_k <- NA_real_
  hill_curve_alpha <- NA_real_
  if (k_tail >= 5 && min(x_tail) > 0) {
    x_sorted <- sort(x_tail[x_tail > 0], decreasing = TRUE)
    k_max <- min(max(10L, floor(length(x_sorted) * tail_fraction)), length(x_sorted) - 1L)
    if (k_max >= 5L) {
      k_seq <- 5L:k_max
      dt_hill <- data.table(log_x = log(x_sorted))
      dt_hill[, csum := cumsum(log_x)]
      denom_vec <- (dt_hill$csum[k_seq] / k_seq) - dt_hill$log_x[k_seq + 1L]
      denom_vec[!is.finite(denom_vec) | denom_vec <= 0] <- NA_real_
      gamma_hat_vec <- denom_vec
      alpha_vec <- 1 / gamma_hat_vec
    } else {
      k_seq <- integer(0)
      alpha_vec <- numeric(0)
    }
    hill_curve_k <- k_seq
    hill_curve_alpha <- alpha_vec
    # default point
    good_idx <- which(is.finite(alpha_vec))
    hill_alpha <- if (length(good_idx) > 0) stats::median(alpha_vec[good_idx], na.rm = TRUE) else NA_real_
  }

  # Power law/Pareto fit via log-log (right tail)
  alpha_pl <- NA_real_
  alpha_pl_se <- NA_real_
  alpha_pl_r2 <- NA_real_
  alpha_pl_p <- NA_real_
  if (k_tail >= 20 && min(x_tail) > 0) {
    # Estimate survival function on tail grid
    xs <- sort(unique(x_tail))
    ec <- stats::ecdf(x_num)
    surv <- 1 - ec(xs)
    keep <- which(surv > 0)
    xs <- xs[keep]
    surv <- surv[keep]
    if (length(xs) >= 20) {
      df_fit <- data.table(x = xs, y = surv)
      df_fit <- df_fit[y <= (tail_fraction + 0.05)]
      if (nrow(df_fit) >= 20) {
        lx <- log(df_fit$x)
        ly <- log(df_fit$y)
        fit <- stats::lm(ly ~ lx)
        alpha_pl <- -as.numeric(stats::coef(fit)["lx"]) # slope with minus sign
        alpha_pl_se <- as.numeric(summary(fit)$coefficients["lx", "Std. Error"])
        alpha_pl_r2 <- as.numeric(summary(fit)$r.squared)
        alpha_pl_p <- as.numeric(summary(fit)$coefficients["lx", "Pr(>|t|)"])
      }
    }
  }

  # Tail severity classification
  tail_label <- "Undefined"
  if (is.finite(kurt_ex) && kurt_ex > 0) tail_label <- "Fat tails (excess kurtosis > 0)"
  if (is.finite(hill_alpha)) {
    if (hill_alpha <= 2) tail_label <- paste0(tail_label, "; infinite variance (alpha<=2)")
    else if (hill_alpha <= 4) tail_label <- paste0(tail_label, "; heavy tails (alpha<=4)")
  }

  # Stats table
  stats_dt <- data.table(
    metric = c("n", "mean", "sd", "var", "median", "mad", "q01", "q05", "q10", "q50", "q90", "q95", "q99",
               "skew", "kurt_excess", "JB_p", "Shapiro_p", "AD_p", "Hill_alpha", "PL_alpha", "PL_alpha_r2", "PL_alpha_p", "tail_label"),
    value = c(n, mean_x, sd_x, var_x, med_x, mad_x, qtls[1], qtls[2], qtls[3], qtls[4], qtls[5], qtls[6], qtls[7],
              skew_ex, kurt_ex, jb_p, shap_p, ad_p, hill_alpha, alpha_pl, alpha_pl_r2, alpha_pl_p, tail_label)
  )

  # Plots
  dens <- stats::density(x_num, n = 2048)
  xs_d <- dens$x
  ys_d <- dens$y
  # normal density with same mean and sd
  ys_n <- stats::dnorm(xs_d, mean = mean_x, sd = sd_x)

  # Pre-calculate bins to avoid passing raw data to plot
  rng <- range(x_num)
  breaks <- seq(rng[1], rng[2], length.out = nbins + 1L)
  binwidth <- (rng[2] - rng[1]) / nbins
  idx <- findInterval(x_num, vec = breaks, rightmost.closed = TRUE, all.inside = TRUE)
  counts <- tabulate(idx, nbins)
  mids <- (breaks[-1L] + breaks[-length(breaks)]) / 2

  # Scale densities to count scale
  scale_factor <- length(x_num) * binwidth
  y_d_count <- ys_d * scale_factor
  y_n_count <- ys_n * scale_factor

  p_hist <- plotly::plot_ly() %>%
    add_bars(x = mids, y = counts, name = "Histogram",
             marker = list(color = "rgba(31,119,180,0.35)"), opacity = 0.6, width = binwidth, hoverinfo = "x+y") %>%
    add_trace(type = "scatter", mode = "lines", x = xs_d, y = y_d_count,
              name = "Kernel density", line = list(color = "rgba(214,39,40,0.8)", width = 2)) %>%
    add_trace(type = "scatter", mode = "lines", x = xs_d, y = y_n_count,
              name = "Normal(mean,sd)", line = list(color = "rgba(44,160,44,0.8)", width = 2, dash = "dash")) %>%
    layout(title = list(text = "Histogram + Densities"), barmode = "overlay")

  ec <- stats::ecdf(x_num)
  xs_ec <- sort(unique(x_num))
  ys_ec <- ec(xs_ec)
  p_ecdf <- plotly::plot_ly(type = "scatter", mode = "lines", x = xs_ec, y = ys_ec, name = "ECDF",
                            line = list(color = "rgba(148,103,189,0.9)", width = 2)) %>%
    layout(title = list(text = "Empirical CDF"), yaxis = list(range = c(0, 1)))

  # QQ to normal (y=x line as separate trace without points)
  q_probs <- seq(0.01, 0.99, by = 0.01)
  q_emp <- stats::quantile(x_num, probs = q_probs, names = FALSE)
  q_th <- stats::qnorm(q_probs, mean = mean_x, sd = sd_x)
  p_qq <- plotly::plot_ly() %>%
    add_trace(type = "scatter", mode = "lines", x = q_th, y = q_th, name = "y=x",
              line = list(color = "rgba(99,99,99,0.6)", dash = "dash")) %>%
    add_trace(type = "scatter", mode = "markers", x = q_th, y = q_emp, name = "QQ",
              marker = list(color = "rgba(255,127,14,0.7)")) %>%
    layout(title = list(text = "QQ to Normal(mean,sd)"), xaxis = list(title = "Theoretical quantiles"), yaxis = list(title = "Empirical quantiles"))

  # Log-log tail and power law fit (reference line without points)
  p_tail <- NULL
  if (k_tail >= 20 && min(x_tail) > 0) {
    xs <- sort(unique(x_tail))
    surv <- 1 - ec(xs)
    keep <- which(surv > 0)
    xs <- xs[keep]
    surv <- surv[keep]
    p_tail <- plotly::plot_ly()
    if (is.finite(alpha_pl)) {
      # Line: log S(x) = a - alpha*log x
      a_hat <- as.numeric(mean(log(surv) + alpha_pl * log(xs)))
      ys_hat <- a_hat - alpha_pl * log(xs)
      p_tail <- p_tail %>% add_trace(type = "scatter", mode = "lines", x = log(xs), y = ys_hat,
                                     name = paste0("PL alpha=", round(alpha_pl, 2)),
                                     line = list(color = "rgba(214,39,40,0.9)", width = 2))
    }
    p_tail <- p_tail %>% add_trace(type = "scatter", mode = "markers", x = log(xs), y = log(surv),
                                   name = "log-log tail", marker = list(color = "rgba(31,119,180,0.7)", size = 6)) %>%
      layout(title = list(text = "Log-Log Tail (Survival)"), xaxis = list(title = "log(x)"), yaxis = list(title = "log(1-F(x))"))
  } else {
    p_tail <- plotly::plot_ly(type = "scatter", mode = "text", x = 0, y = 0, text = "Insufficient data for tail",
                              textposition = "middle center") %>% layout(title = list(text = "Log-Log Tail"))
  }

  # Hill plot (downsampling to <= 200 points)
  p_hill <- NULL
  if (is.numeric(hill_curve_k) && length(hill_curve_k) > 1 && any(is.finite(hill_curve_alpha))) {
    k_plot <- hill_curve_k
    a_plot <- hill_curve_alpha
    if (length(k_plot) > 200L) {
      idx <- unique(as.integer(round(seq(1, length(k_plot), length.out = 200L))))
      k_plot <- k_plot[idx]
      a_plot <- a_plot[idx]
    }
    p_hill <- plotly::plot_ly(type = "scatter", mode = "lines+markers", x = k_plot, y = a_plot,
                              name = "alpha(k)", line = list(color = "rgba(23,190,207,0.8)"), marker = list(size = 6)) %>%
      layout(title = list(text = "Hill alpha(k)"), xaxis = list(title = "k (upper order statistics)"), yaxis = list(title = "alpha"))
  } else {
    p_hill <- plotly::plot_ly(type = "scatter", mode = "text", x = 0, y = 0, text = "Hill unavailable",
                              textposition = "middle center") %>% layout(title = list(text = "Hill"))
  }

  # Violin/Box without raw data: density "violin" and quantile/whisker lines
  q2575 <- stats::quantile(x_num, probs = c(0.25, 0.5, 0.75), names = FALSE)
  q25 <- q2575[1]; q50 <- q2575[2]; q75 <- q2575[3]
  iqr <- q75 - q25
  whisk_lo <- max(min(x_num), q25 - 1.5 * iqr)
  whisk_hi <- min(max(x_num), q75 + 1.5 * iqr)
  w_max <- 0.5
  ymax <- max(ys_d)
  ys_norm <- if (is.finite(ymax) && ymax > 0) ys_d / ymax else ys_d
  w <- ys_norm * w_max
  x_poly <- c(-w, rev(w))
  y_poly <- c(xs_d, rev(xs_d))
  p_box <- plotly::plot_ly() %>%
    add_trace(type = "scatter", mode = "lines", x = x_poly, y = y_poly, fill = "toself", name = "Violin",
              line = list(color = "rgba(148,103,189,1)"), fillcolor = "rgba(148,103,189,0.35)", opacity = 0.6, hoverinfo = "skip") %>%
    add_segments(x = -w_max, xend = w_max, y = q25, yend = q25, name = "Q1",
                 line = list(color = "rgba(99,99,99,0.7)", width = 2)) %>%
    add_segments(x = -w_max, xend = w_max, y = q50, yend = q50, name = "Median",
                 line = list(color = "rgba(214,39,40,0.9)", width = 2)) %>%
    add_segments(x = -w_max, xend = w_max, y = q75, yend = q75, name = "Q3",
                 line = list(color = "rgba(99,99,99,0.7)", width = 2)) %>%
    add_segments(x = 0, xend = 0, y = whisk_lo, yend = q25, name = "Whisker-",
                 line = list(color = "rgba(99,99,99,0.5)", width = 2, dash = "dot")) %>%
    add_segments(x = 0, xend = 0, y = q75, yend = whisk_hi, name = "Whisker+",
                 line = list(color = "rgba(99,99,99,0.5)", width = 2, dash = "dot")) %>%
    layout(title = list(text = "Violin/Box"), xaxis = list(showticklabels = FALSE, zeroline = FALSE, title = ""))

  # # Stats table
  # tbl <- plotly::plot_ly(type = "table",
  #                        header = list(values = c("Metric", "Value"), align = c("left", "right")),
  #                        cells = list(values = list(stats_dt$metric, as.character(stats_dt$value)), align = c("left", "right")))

  fig <- subplot(list(p_hist, p_ecdf, p_qq, p_tail, p_hill, p_box), nrows = 3, shareX = FALSE, shareY = FALSE, margin = 0.05, titleX = TRUE, titleY = TRUE)
  fig <- fig %>% layout(title = list(text = tail_label))

  result <- list(plot = fig, stats = stats_dt, tests = list(Shapiro_p = shap_p, AD_p = ad_p, JB_p = jb_p),
                 tail = list(Hill_alpha = hill_alpha, PL_alpha = alpha_pl, PL_r2 = alpha_pl_r2, PL_p = alpha_pl_p, label = tail_label))

  if (!is.null(save_qs_path)) {
    try(qs::qsave(result, file = save_qs_path), silent = TRUE)
  }
  # mean, sd, var: mean/standard deviation/variance. Unreliable with fat tails (may "float"). Compare with median/MAD.
  # median, MAD: robust center/scale. If significantly different from mean/sd - tails/outliers are substantial.
  # quantiles (q01...q99): basis for risk; compare right/left tails, dynamics across windows.
  # skew: asymmetry. >0 - longer right tail; <0 - left tail.
  # kurt_excess: excess kurtosis (0 for normal). 0.2-1: moderately fat; 1-3: heavy; >3: extreme tails.
  # JB_p: Jarque-Bera. <0.05 - not normal (sensitive to skewness/excess kurtosis).
  # Shapiro_p: <0.05 - not normal; with large n almost always rejects - look at tail metrics.
  # AD_p: Anderson-Darling, more "tail-sensitive". <0.05 - not normal.
  # Hill_alpha: tail index α (right). α≤2 - infinite variance; 2<α≤4 - heavy tails; α>4 - closer to "thin". Evaluate not a point but plateau on Hill plot.
  # PL_alpha: power law exponent from log-survival regression; α estimate. Trust with good linear segment and high consistency with Hill.
  # PL_alpha_r2: quality of power law fit on tail. >0.9 - acceptable; lower - interpret with caution.
  # PL_alpha_p: slope significance; secondary to R² and stability.
  # tail_label: brief verbal classification of tails based on combination of indicators.
  return(result)
}

#### Usage Examples ####

# Normal distribution
normal_example <- BuildDistDashboard(rnorm(10000))

# Strange distribution with local bumps (mixture of normals + noise)
strange_example <- BuildDistDashboard(
  c(rnorm(3000, mean = 0, sd = 1),
    rnorm(2000, mean = 3, sd = 0.5),
    rnorm(1500, mean = -2, sd = 0.8),
    rnorm(1000, mean = 6, sd = 0.3),
    rnorm(500, mean = -4, sd = 0.4),
    rnorm(1000, mean = 0, sd = 0.1) + runif(1000, -0.5, 0.5))
)

# Pareto distribution (heavy tails)
pareto_example <- BuildDistDashboard(Pareto::rPareto(10000, alpha = 1.5, t = 0.05))

# Mixed normal + Pareto (realistic financial returns)
mixed_example <- BuildDistDashboard(rnorm(10000) + runif(10000, -1, 1) * Pareto::rPareto(10000, alpha = 1.5, t = 0.05))

# Cauchy distribution (very heavy tails)
cauchy_example <- BuildDistDashboard(rcauchy(10000, location = 0, scale = 1))

# Levy distribution (stable, alpha=0.5)
levy_example <- BuildDistDashboard(stabledist::rstable(10000, alpha = 0.5, beta = 0, gamma = 1, delta = 0))
