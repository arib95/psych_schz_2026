suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

DX_FILE <- "data/wide_diagnoses.csv"
PSY_FILE <- "data/psychometric_matrix.csv"
OUTCOME_VAR <- "SCID.DIAG.Schizophrenia"
ID_VAR <- "participant_id"
OUT_DIR <- "out"

K_FOLDS <- 5L
SEED <- 42L
BOOT_B <- 1000L
BOOT_SEED <- SEED + 1000L

N_PERM <- 100L

DOM_ENABLED <- TRUE
DOM_TOP_K <- 19L
DOM_MAXIT <- 50L
DOM_CORES <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

PREDICTORS <- c(
  "psy.tci.tci82t",
  "psy.tci.tci34t",
  "psy.tci.tci203t",
  "psy.tci.tci148t",
  "psy.tci.tci141t",
  "psy.mpq.mpq47",
  "psy.mpq.mpq209",
  "psy.mpq.mpq2",
  "psy.mpq.mpq172",
  "psy.mpq.mpq147",
  "psy.eysenck.eysenck44",
  "psy.eysenck.eysenck19",
  "psy.dickman.dick40",
  "psy.chapper.chapper22",
  "psy.chapper.chapper14",
  "psy.chapinf.chapinf12",
  
  "psy.tci.tci109t",
  "psy.mpq.mpq92",
  "psy.mpq.mpq219"
)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

msg <- function(fmt, ...) cat(sprintf(paste0(fmt, "\n"), ...))

`%||%` <- function(x, y) if (is.null(x)) y else x

binom_wilson_ci <- function(x, n, conf = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0) return(c(lo = NA_real_, hi = NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  den <- 1 + z^2 / n
  ctr <- (phat + z^2 / (2 * n)) / den
  rad <- (z * sqrt((phat * (1 - phat) / n) + z^2 / (4 * n^2))) / den
  c(lo = max(0, ctr - rad), hi = min(1, ctr + rad))
}

char_to_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "na", "N/A", ".")] <- NA_character_
  suppressWarnings(out <- as.numeric(x))
  if (all(is.na(out)) && any(!is.na(x))) {
    xl <- tolower(x)
    out <- rep(NA_real_, length(xl))
    out[xl %in% c("yes", "y", "true", "t", "case", "present", "1")] <- 1
    out[xl %in% c("no", "n", "false", "f", "control", "absent", "0")] <- 0
  }
  out
}

coerce_binary01 <- function(x) {
  x_num <- char_to_num(x)
  ux <- sort(unique(stats::na.omit(x_num)))
  if (length(ux) == 0L) return(x_num)
  if (all(ux %in% c(0, 1))) return(x_num)
  if (length(ux) == 2L) {
    out <- rep(NA_real_, length(x_num))
    out[x_num == ux[1]] <- 0
    out[x_num == ux[2]] <- 1
    return(out)
  }
  stop(sprintf("Outcome '%s' is not binary after coercion.", OUTCOME_VAR))
}

safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

calc_auc <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- as.numeric(y[keep])
  score <- as.numeric(score[keep])
  if (length(y) < 2L || length(unique(y)) < 2L) return(NA_real_)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

logloss <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}


brier_point <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- as.numeric(y[keep])
  p <- as.numeric(p[keep])
  if (!length(y)) return(NA_real_)
  mean((y - p)^2)
}

calibration_metrics_binary <- function(y, p, min_bin = 25L, max_bins = 8L) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  
  if (length(unique(y)) < 2L || length(y) < 30L) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      brier = NA_real_,
      ece = NA_real_,
      n_bins = 0L,
      n_total = length(y),
      curve = data.frame(
        p_hat = numeric(0),
        obs   = numeric(0),
        n     = integer(0)
      ),
      line = data.frame(
        p_hat = numeric(0),
        p     = numeric(0),
        lo    = numeric(0),
        hi    = numeric(0)
      )
    ))
  }
  
  p_clip <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  lp <- qlogis(p_clip)
  
  fit_cal <- try(glm(y ~ lp, family = binomial), silent = TRUE)
  
  intercept <- NA_real_
  slope <- NA_real_
  line <- data.frame(
    p_hat = numeric(0),
    p     = numeric(0),
    lo    = numeric(0),
    hi    = numeric(0)
  )
  
  if (!inherits(fit_cal, "try-error")) {
    cf <- coef(fit_cal)
    if (length(cf) >= 2L && all(is.finite(cf[1:2]))) {
      intercept <- unname(cf[1])
      slope <- unname(cf[2])
    }
    
    gp <- seq(0, 1, length.out = 200L)
    gp_clip <- pmin(pmax(gp, 1e-6), 1 - 1e-6)
    
    pr <- predict(
      fit_cal,
      newdata = data.frame(lp = qlogis(gp_clip)),
      type = "link",
      se.fit = TRUE
    )
    
    line <- data.frame(
      p_hat = gp,
      p  = plogis(pr$fit),
      lo = plogis(pr$fit - 1.96 * pr$se.fit),
      hi = plogis(pr$fit + 1.96 * pr$se.fit)
    )
  }
  
  brier <- mean((p_clip - y)^2)
  
  nb <- max(3L, min(as.integer(max_bins), floor(length(p_clip) / min_bin)))
  br <- unique(as.numeric(stats::quantile(
    p_clip,
    probs = seq(0, 1, length.out = nb + 1L),
    na.rm = TRUE,
    names = FALSE
  )))
  if (length(br) < 3L) br <- c(0, 0.5, 1)
  
  g <- cut(p_clip, breaks = br, include.lowest = TRUE)
  
  T <- data.frame(
    p_hat = as.numeric(tapply(p_clip, g, mean)),
    obs   = as.numeric(tapply(y,      g, mean)),
    n     = as.integer(tapply(y,      g, length))
  )
  
  T <- T[is.finite(T$p_hat) & is.finite(T$obs) & is.finite(T$n), , drop = FALSE]
  
  list(
    intercept = intercept,
    slope     = slope,
    brier     = brier,
    ece       = weighted.mean(abs(T$p_hat - T$obs), w = T$n),
    n_bins    = nrow(T),
    n_total   = length(y),
    curve     = T,
    line      = line
  )
}

compute_binary_metrics_with_ci <- function(y01, p_hat, boot_B = BOOT_B, boot_seed = BOOT_SEED) {
  y01 <- as.integer(y01 > 0)
  p_hat <- pmin(pmax(as.numeric(p_hat), 1e-6), 1 - 1e-6)
  keep <- is.finite(y01) & is.finite(p_hat)
  y01 <- y01[keep]
  p_hat <- p_hat[keep]
  n <- length(y01)
  if (n < 5L || length(unique(y01)) < 2L) {
    return(list(
      auc = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_,
      auprc = NA_real_, auprc_lo = NA_real_, auprc_hi = NA_real_,
      auprg = NA_real_, auprg_lo = NA_real_, auprg_hi = NA_real_,
      logloss = NA_real_, logloss_lo = NA_real_, logloss_hi = NA_real_,
      brier = NA_real_, brier_lo = NA_real_, brier_hi = NA_real_,
      cal_int = NA_real_, cal_int_lo = NA_real_, cal_int_hi = NA_real_,
      cal_slope = NA_real_, cal_slope_lo = NA_real_, cal_slope_hi = NA_real_,
      cal_ece = NA_real_, cal_ece_lo = NA_real_, cal_ece_hi = NA_real_,
      cal_curve = data.frame(p_hat = numeric(0), obs = numeric(0), n = integer(0), lo = numeric(0), hi = numeric(0)),
      cal_line = data.frame(p_hat = numeric(0), p = numeric(0), lo = numeric(0), hi = numeric(0))
    ))
  }
  calib <- calibration_metrics_binary(y01, p_hat)
  ci_auc <- boot_ci_balanced(y01, p_hat, auc_point, B = boot_B, seed = boot_seed + 1L)
  ci_auprc <- boot_ci_balanced(y01, p_hat, auprc_point, B = boot_B, seed = boot_seed + 2L)
  ci_auprg <- boot_ci_balanced(y01, p_hat, function(yy, pp) prg_flach(yy, pp)$auprg, B = boot_B, seed = boot_seed + 3L)  
  ci_logloss <- boot_ci_balanced(y01, p_hat, logloss_binary_point, B = boot_B, seed = boot_seed + 4L)
  ci_brier <- boot_ci_balanced(y01, p_hat, brier_point, B = boot_B, seed = boot_seed + 5L)
  ci_cal_int <- boot_ci_balanced(y01, p_hat, function(yy, pp) calibration_metrics_binary(yy, pp)$intercept, B = max(200L, boot_B %/% 2L), seed = boot_seed + 6L)
  ci_cal_slope <- boot_ci_balanced(y01, p_hat, function(yy, pp) calibration_metrics_binary(yy, pp)$slope, B = max(200L, boot_B %/% 2L), seed = boot_seed + 7L)
  ci_ece <- boot_ci_balanced(y01, p_hat, function(yy, pp) calibration_metrics_binary(yy, pp)$ece, B = max(200L, boot_B %/% 2L), seed = boot_seed + 8L)
  list(
    auc = unname(ci_auc[['point']]), auc_lo = unname(ci_auc[['lo']]), auc_hi = unname(ci_auc[['hi']]),
    auprc = unname(ci_auprc[['point']]), auprc_lo = unname(ci_auprc[['lo']]), auprc_hi = unname(ci_auprc[['hi']]),
    auprg = unname(ci_auprg[['point']]), auprg_lo = unname(ci_auprg[['lo']]), auprg_hi = unname(ci_auprg[['hi']]),
    logloss = unname(ci_logloss[['point']]), logloss_lo = unname(ci_logloss[['lo']]), logloss_hi = unname(ci_logloss[['hi']]),
    brier = unname(ci_brier[['point']]), brier_lo = unname(ci_brier[['lo']]), brier_hi = unname(ci_brier[['hi']]),
    cal_int = calib$intercept, cal_int_lo = unname(ci_cal_int[['lo']]), cal_int_hi = unname(ci_cal_int[['hi']]),
    cal_slope = calib$slope, cal_slope_lo = unname(ci_cal_slope[['lo']]), cal_slope_hi = unname(ci_cal_slope[['hi']]),
    cal_ece = calib$ece, cal_ece_lo = unname(ci_ece[['lo']]), cal_ece_hi = unname(ci_ece[['hi']]),
    cal_curve = calib$curve,
    cal_line = calib$line
  )
}

boot_ci_balanced <- function(y, p, FUN, B = 1000L, seed = 42L) {
  set.seed(seed)
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  i0 <- which(y == 0L)
  i1 <- which(y == 1L)
  pt <- suppressWarnings(as.numeric(FUN(y, p)))
  
  if (!length(i0) || !length(i1) || !is.finite(pt)) {
    return(c(point = pt, lo = NA, hi = NA))
  }
  
  vals <- replicate(B, {
    ii <- c(sample(i0, length(i0), TRUE), sample(i1, length(i1), TRUE))
    suppressWarnings(as.numeric(FUN(y[ii], p[ii])))
  })
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(c(point = pt, lo = NA, hi = NA))
  }
  
  c(point = pt, lo = as.numeric(quantile(vals, 0.025)), hi = as.numeric(quantile(vals, 0.975)))
}

auc_point <- calc_auc
logloss_binary_point <- logloss

roc_curve_df <- function(y, p) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  
  P <- sum(y == 1L)
  N <- sum(y == 0L)
  if (length(y) < 2L || P == 0L || N == 0L) {
    return(data.frame(
      threshold = numeric(0),
      fpr = numeric(0),
      tpr = numeric(0),
      specificity = numeric(0),
      sensitivity = numeric(0)
    ))
  }
  
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  p <- p[o]
  
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  
  tpr <- tp / P
  fpr <- fp / N
  
  data.frame(
    threshold   = c(Inf, p, -Inf),
    fpr         = c(0, fpr, 1),
    tpr         = c(0, tpr, 1),
    specificity = c(1, 1 - fpr, 0),
    sensitivity = c(0, tpr, 1)
  )
}

pr_curve_df <- function(y, p) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  
  P <- sum(y == 1L)
  N <- sum(y == 0L)
  if (length(y) < 2L || P == 0L || N == 0L) {
    return(data.frame(
      threshold = numeric(0),
      recall = numeric(0),
      precision = numeric(0)
    ))
  }
  
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  p <- p[o]
  
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  
  recall <- tp / P
  precision <- tp / pmax(tp + fp, 1L)
  
  data.frame(
    threshold = c(Inf, p),
    recall    = c(0, recall),
    precision = c(mean(y), precision)
  )
}

auprc_point <- function(y, p) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  n <- length(y)
  P <- sum(y == 1L)
  N <- n - P
  if (n < 2L || P == 0L || N == 0L) {
    return(NA_real_)
  }
  
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  prec <- tp / pmax(tp + fp, 1L)
  
  sum(prec[y == 1L]) / P
}

prg_flach <- function(y, p, eps = 1e-12) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(y) == 0L || length(unique(y)) < 2L) {
    return(list(auprg = NA_real_))
  }
  
  pi0 <- mean(y == 1L)
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  n1 <- sum(y == 1L)
  
  R  <- tp / pmax(n1, 1L)
  Pp <- tp / pmax(tp + fp, 1L)
  
  R  <- c(0, R, 1)
  Pp <- c(pi0, Pp, pi0)
  
  PR <- data.frame(R = R, P = Pp)
  PR <- PR[order(PR$R, -PR$P), , drop = FALSE]
  PR <- aggregate(P ~ R, PR, max)
  PR <- PR[order(PR$R), , drop = FALSE]
  
  recG  <- (PR$R - pi0) / ((1 - pi0) * pmax(PR$R, eps))
  precG <- (PR$P - pi0) / ((1 - pi0) * pmax(PR$P, eps))
  
  recG  <- pmin(pmax(recG,  0), 1)
  precG <- pmin(pmax(precG, 0), 1)
  
  o2 <- order(recG)
  recG  <- recG[o2]
  precG <- precG[o2]
  
  auprg <- if (length(recG) >= 2L) {
    sum(0.5 * (precG[-1] + precG[-length(precG)]) * diff(recG))
  } else {
    NA_real_
  }
  
  list(auprg = as.numeric(auprg))
}

prg_curve_df <- function(y, p, eps = 1e-12) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  
  if (length(y) < 2L || length(unique(y)) < 2L) {
    return(data.frame(
      threshold = numeric(0),
      recall_gain = numeric(0),
      precision_gain = numeric(0)
    ))
  }
  
  pi0 <- mean(y == 1L)
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  p <- p[o]
  
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  P <- sum(y == 1L)
  
  recall <- tp / pmax(P, 1L)
  precision <- tp / pmax(tp + fp, 1L)
  
  rec_gain  <- (recall    - pi0) / ((1 - pi0) * pmax(recall,    eps))
  prec_gain <- (precision - pi0) / ((1 - pi0) * pmax(precision, eps))
  
  out <- data.frame(
    threshold = p,
    recall_gain = rec_gain,
    precision_gain = prec_gain
  )
  
  out$recall_gain[!is.finite(out$recall_gain)] <- 0
  out$precision_gain[!is.finite(out$precision_gain)] <- 0
  
  out$recall_gain <- pmin(pmax(out$recall_gain, 0), 1)
  out$precision_gain <- pmin(pmax(out$precision_gain, 0), 1)
  
  out <- rbind(
    data.frame(threshold = Inf,  recall_gain = 0, precision_gain = 1),
    out,
    data.frame(threshold = -Inf, recall_gain = 1, precision_gain = 0)
  )
  
  out <- out[order(out$recall_gain, -out$precision_gain), , drop = FALSE]
  out <- aggregate(precision_gain ~ recall_gain, data = out, FUN = max)
  out <- out[order(out$recall_gain), , drop = FALSE]
  out
}

auprg_point <- function(y, p) {
  prg_flach(y, p)$auprg
}

save_plot_gg_local <- function(plot, filename, width, height, dpi = 300,
                               save_pdf = TRUE, save_rds = FALSE) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = "white"
  )
  
  if (isTRUE(save_pdf)) {
    pdf_file <- sub("\\.[A-Za-z0-9]+$", ".pdf", filename)
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      device = grDevices::cairo_pdf,
      units = "in",
      bg = "white"
    )
  }
  
  if (isTRUE(save_rds)) {
    rds_file <- sub("\\.[A-Za-z0-9]+$", ".rds", filename)
    saveRDS(plot, file = rds_file)
  }
  
  invisible(plot)
}

theme_pub_local <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(size = base_size),
      legend.text      = ggplot2::element_text(size = base_size - 1)
    )
}

write_binary_curves <- function(y01, p_hat, out_png,
                                title = "",
                                subtitle = "",
                                metrics = NULL) {
  y01   <- as.integer(y01 > 0)
  p_hat <- as.numeric(p_hat)
  ok <- is.finite(y01) & is.finite(p_hat)
  y01   <- y01[ok]
  p_hat <- pmin(pmax(p_hat[ok], 1e-6), 1 - 1e-6)
  
  if (length(y01) < 2L || length(unique(y01)) < 2L) {
    msg("[curves] skipped %s: need both outcome classes.", out_png)
    return(invisible(NULL))
  }
  
  if (is.null(metrics)) {
    metrics <- compute_binary_metrics_with_ci(y01, p_hat)
  }
  
  df_roc <- roc_curve_df(y01, p_hat)
  df_pr  <- pr_curve_df(y01, p_hat)
  df_prg <- prg_curve_df(y01, p_hat)
  prev_v <- mean(y01 == 1L)
  
  theme_c <- theme_pub_local(11)
  
  auc_txt <- if (is.finite(metrics$auc_lo) && is.finite(metrics$auc_hi)) {
    sprintf("AUC = %.3f [%.3f, %.3f]", metrics$auc, metrics$auc_lo, metrics$auc_hi)
  } else {
    sprintf("AUC = %.3f", metrics$auc)
  }
  
  auprc_txt <- if (is.finite(metrics$auprc_lo) && is.finite(metrics$auprc_hi)) {
    sprintf("AUPRC = %.3f [%.3f, %.3f]", metrics$auprc, metrics$auprc_lo, metrics$auprc_hi)
  } else {
    sprintf("AUPRC = %.3f", metrics$auprc)
  }
  
  auprg_txt <- if (is.finite(metrics$auprg_lo) && is.finite(metrics$auprg_hi)) {
    sprintf("AUPRG = %.3f [%.3f, %.3f]", metrics$auprg, metrics$auprg_lo, metrics$auprg_hi)
  } else {
    sprintf("AUPRG = %.3f", metrics$auprg)
  }
  
  p_roc <- ggplot2::ggplot(df_roc, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey70", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(
      title = paste("ROC |", auc_txt),
      x = "False positive rate",
      y = "True positive rate"
    ) +
    theme_c
  
  p_pr <- ggplot2::ggplot(df_pr, ggplot2::aes(x = recall, y = precision)) +
    ggplot2::geom_hline(yintercept = prev_v, linetype = 2, colour = "grey70", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(
      title = paste("PRC |", auprc_txt),
      x = "Recall",
      y = "Precision"
    ) +
    theme_c
  
  p_prg <- if (nrow(df_prg) && is.finite(metrics$auprg)) {
    ggplot2::ggplot(df_prg, ggplot2::aes(x = recall_gain, y = precision_gain)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey70", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey70", linewidth = 0.5) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::labs(
        title = paste("PRG |", auprg_txt),
        x = "Recall gain",
        y = "Precision gain"
      ) +
      theme_c
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::labs(title = "PRG unavailable", x = NULL, y = NULL) +
      theme_c
  }
  
  ttl <- paste(c(title, subtitle), collapse = "\n")
  ttl <- gsub("\n+$", "", ttl)
  
  p_all <- (p_roc | p_pr | p_prg) +
    patchwork::plot_annotation(title = ttl)
  
  save_plot_gg_local(p_all, out_png, width = 12.5, height = 3.85, dpi = 300)
  
  cal_png <- sub("\\.png$", "_calibration.png", out_png, ignore.case = TRUE)
  
  if (!is.null(metrics$cal_curve) && nrow(metrics$cal_curve) > 0) {
    df_cal  <- metrics$cal_curve
    line_df <- metrics$cal_line
    
    cal_txt <- sprintf(
      "Calib: intercept = %.3f | slope = %.3f | Brier = %.3f | ECE = %.3f (bins = %d)",
      metrics$cal_int, metrics$cal_slope, metrics$brier, metrics$cal_ece, nrow(df_cal)
    )
    
    cap_txt <- NULL
    if (isTRUE(nrow(df_cal) > 0) && isTRUE(sum(df_cal$n) > 0)) {
      cap_txt <- sprintf(
        "%d bins of ~%.1f patients each (n=%d)",
        nrow(df_cal), sum(df_cal$n) / nrow(df_cal), sum(df_cal$n)
      )
    }
    
    df_cal$p_hat <- pmin(pmax(df_cal$p_hat, 0), 1)
    df_cal$obs   <- pmin(pmax(df_cal$obs,   0), 1)
    
    p_cal <- ggplot2::ggplot(df_cal, ggplot2::aes(x = p_hat, y = obs)) +
      ggplot2::geom_abline(
        intercept = 0, slope = 1,
        linetype = 2, colour = "grey70", linewidth = 0.5
      )
    
    if (!is.null(line_df) && nrow(line_df) > 0) {
      p_cal <- p_cal +
        ggplot2::geom_ribbon(
          data = line_df,
          ggplot2::aes(x = p_hat, ymin = lo, ymax = hi),
          inherit.aes = FALSE,
          alpha = 0.15
        ) +
        ggplot2::geom_line(
          data = line_df,
          ggplot2::aes(x = p_hat, y = p),
          inherit.aes = FALSE,
          linewidth = 0.9
        )
    }
    
    p_cal <- p_cal +
      ggplot2::geom_point(shape = 16, size = 1.8, alpha = 0.7) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::labs(
        x = "Predicted probability",
        y = "Observed event rate",
        title = paste("Calibration —", title),
        subtitle = cal_txt,
        caption = cap_txt
      ) +
      theme_pub_local(11) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "none"
      )
    
    save_plot_gg_local(p_cal, cal_png, width = 5.5, height = 5, dpi = 300)
  }
  
  invisible(metrics)
}

logit_or <- function(y, x) {
  keep <- is.finite(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]
  if (length(y) < 10L || length(unique(y)) < 2L || stats::sd(x) <= 0) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, beta = NA_real_))
  }
  fit <- try(stats::glm(y ~ x, family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, beta = NA_real_))
  }
  cf <- summary(fit)$coefficients
  beta <- cf["x", "Estimate"]
  se   <- cf["x", "Std. Error"]
  p    <- cf["x", "Pr(>|z|)"]
  list(or = exp(beta), lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se), p = p, beta = beta)
}

make_folds_bin <- function(y, K = 5L, seed = 1L) {
  set.seed(seed)
  y <- as.integer(y > 0)
  fid <- integer(length(y))
  i1 <- which(y == 1L)
  i0 <- which(y == 0L)
  fid[i1] <- sample(rep(seq_len(K), length.out = length(i1)))
  fid[i0] <- sample(rep(seq_len(K), length.out = length(i0)))
  fid
}

oof_ridge_bin <- function(X, y, K = 5L, seed = 1L, n_perm = 0L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for this script.")
  }
  
  y <- as.integer(y > 0)
  X <- as.matrix(X)
  fid <- make_folds_bin(y, K = K, seed = seed)
  
  p_oof <- rep(NA_real_, length(y))
  eta_oof <- rep(NA_real_, length(y))
  coef_list <- list()
  perm_list <- list()
  
  for (k in sort(unique(fid))) {
    msg("[ridge] outer fold %d / %d", k, K)
    
    tr <- fid != k
    te <- fid == k
    
    Xtr <- X[tr, , drop = FALSE]
    Xte <- X[te, , drop = FALSE]
    ytr <- y[tr]
    yte <- y[te]
    
    mu <- colMeans(Xtr, na.rm = TRUE)
    sdv <- apply(Xtr, 2, stats::sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    
    Xtrs <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
    Xtes <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")
    
    foldid_inner <- make_folds_bin(ytr, K = 5L, seed = seed + k)
    fit <- glmnet::cv.glmnet(
      x = Xtrs,
      y = ytr,
      alpha = 0,
      family = "binomial",
      foldid = foldid_inner,
      parallel = FALSE
    )
    
    p_te <- as.numeric(stats::predict(fit, newx = Xtes, s = "lambda.min", type = "response"))
    eta_te <- as.numeric(stats::predict(fit, newx = Xtes, s = "lambda.min", type = "link"))
    
    p_oof[te] <- p_te
    eta_oof[te] <- eta_te
    
    cf <- as.matrix(stats::coef(fit, s = "lambda.min"))[, 1]
    nm <- rownames(as.matrix(stats::coef(fit, s = "lambda.min")))
    names(cf) <- nm
    coef_list[[length(coef_list) + 1L]] <- cf
    
    if (n_perm > 0L) {
      base_auc <- calc_auc(yte, p_te)
      base_ll  <- logloss(yte, p_te)
      
      perm_fold <- vector("list", n_perm)
      
      for (b in seq_len(n_perm)) {
        perm_dt <- data.table(
          item = colnames(X),
          dAUC = NA_real_,
          dLogLoss = NA_real_,
          fold = k,
          perm_rep = b
        )
        
        for (j in seq_len(ncol(X))) {
          xperm <- Xtes
          xperm[, j] <- sample(xperm[, j])
          p_perm <- as.numeric(stats::predict(fit, newx = xperm, s = "lambda.min", type = "response"))
          perm_dt[j, `:=`(
            dAUC = base_auc - calc_auc(yte, p_perm),
            dLogLoss = logloss(yte, p_perm) - base_ll
          )]
        }
        
        perm_fold[[b]] <- perm_dt
      }
      
      perm_list[[length(perm_list) + 1L]] <- rbindlist(perm_fold, use.names = TRUE, fill = TRUE)
    }
  }
  
  coef_names <- unique(unlist(lapply(coef_list, names)))
  coef_mat <- matrix(
    NA_real_,
    nrow = length(coef_names),
    ncol = length(coef_list),
    dimnames = list(coef_names, paste0("fold", seq_along(coef_list)))
  )
  
  for (j in seq_along(coef_list)) {
    cj <- coef_list[[j]]
    coef_mat[names(cj), j] <- cj
  }
  
  perm_summary <- NULL
  if (length(perm_list)) {
    perm_all <- rbindlist(perm_list, use.names = TRUE, fill = TRUE)
    perm_summary <- perm_all[, .(
      perm_dAUC_mean = mean(dAUC, na.rm = TRUE),
      perm_dAUC_sd = stats::sd(dAUC, na.rm = TRUE),
      perm_dLogLoss_mean = mean(dLogLoss, na.rm = TRUE),
      perm_dLogLoss_sd = stats::sd(dLogLoss, na.rm = TRUE)
    ), by = item]
  }
  
  list(
    p_oof = p_oof,
    eta_oof = eta_oof,
    fid = fid,
    coef_mean = rowMeans(coef_mat, na.rm = TRUE),
    coef_sd = apply(coef_mat, 1L, stats::sd, na.rm = TRUE),
    coef_sign_stability = rowMeans(
      sign(coef_mat) == sign(rowMeans(coef_mat, na.rm = TRUE)),
      na.rm = TRUE
    ),
    perm_summary = perm_summary
  )
}

dominance_mcfadden_exact <- function(X, y, predictors = colnames(X),
                                     ncores = 1L, maxit = 50L) {
  Xmat <- as.matrix(X[, predictors, drop = FALSE])
  storage.mode(Xmat) <- "double"
  yvec <- as.numeric(y)
  p <- ncol(Xmat)
  
  if (p < 1L) stop("No predictors supplied to dominance_mcfadden_exact().")
  if (p > 30L) stop("Bitmask implementation supports at most 30 predictors.")
  
  masks <- 0:(2^p - 1L)
  
  fit0 <- stats::glm.fit(
    x = matrix(1, nrow = nrow(Xmat), ncol = 1L),
    y = yvec,
    family = stats::binomial(),
    control = stats::glm.control(maxit = maxit, epsilon = 1e-8)
  )
  if (!isTRUE(fit0$converged)) {
    stop("Null logistic model did not converge in dominance_mcfadden_exact().")
  }
  
  p0 <- pmin(pmax(fit0$fitted.values, 1e-12), 1 - 1e-12)
  ll0 <- sum(stats::dbinom(yvec, size = 1, prob = p0, log = TRUE))
  
  subset_r2_one <- function(mask, Xmat, yvec, p, ll0, maxit) {
    if (mask == 0L) return(0)
    
    idx <- as.logical(as.integer(intToBits(as.integer(mask)))[seq_len(p)])
    xx <- cbind(1, Xmat[, idx, drop = FALSE])
    
    fit <- try(
      stats::glm.fit(
        x = xx,
        y = yvec,
        family = stats::binomial(),
        control = stats::glm.control(maxit = maxit, epsilon = 1e-8)
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
      return(NA_real_)
    }
    
    pp <- pmin(pmax(fit$fitted.values, 1e-12), 1 - 1e-12)
    ll1 <- sum(stats::dbinom(yvec, size = 1, prob = pp, log = TRUE))
    
    1 - (ll1 / ll0)
  }
  
  worker_fun <- function(mask) {
    tryCatch(
      subset_r2_one(mask, Xmat = Xmat, yvec = yvec, p = p, ll0 = ll0, maxit = maxit),
      error = function(e) structure(
        NA_real_,
        class = "dom_worker_error",
        worker_message = conditionMessage(e)
      )
    )
  }
  
  msg("[dom] exact McFadden dominance on p=%d predictors (~%d subset models) using %d cores",
      p, length(masks) - 1L, ncores)
  
  t0 <- proc.time()[3]
  
  if (ncores <= 1L) {
    res <- lapply(masks, worker_fun)
  } else {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = c("Xmat", "yvec", "p", "ll0", "maxit", "subset_r2_one", "worker_fun"),
      envir = environment()
    )
    
    res <- parallel::parLapplyLB(cl, masks, worker_fun)
  }
  
  is_err <- vapply(res, inherits, logical(1), what = "dom_worker_error")
  if (any(is_err)) {
    ex <- res[[which(is_err)[1L]]]
    stop(sprintf(
      "[dom] parallel worker error in dominance calculation: %s",
      attr(ex, "worker_message")
    ))
  }
  
  r2 <- as.numeric(unlist(res, use.names = FALSE))
  names(r2) <- as.character(masks)
  
  subset_size <- vapply(
    masks,
    function(m) sum(as.integer(intToBits(as.integer(m)))[seq_len(p)]),
    integer(1)
  )
  
  dom <- numeric(p)
  names(dom) <- predictors
  
  for (j in seq_len(p)) {
    bit <- bitwShiftL(1L, j - 1L)
    
    base_masks <- masks[bitwAnd(masks, bit) == 0L]
    plus_masks <- bitwOr(base_masks, bit)
    
    base_idx <- match(as.character(base_masks), names(r2))
    plus_idx <- match(as.character(plus_masks), names(r2))
    
    deltas <- r2[plus_idx] - r2[base_idx]
    
    lev <- subset_size[match(base_masks, masks)]
    lev_means <- tapply(deltas, lev, function(z) mean(z, na.rm = TRUE))
    dom[j] <- mean(lev_means[is.finite(lev_means)])
  }
  
  elapsed <- (proc.time()[3] - t0) / 60
  usable <- sum(is.finite(r2))
  
  msg("[dom] finished in %.2f minutes; %d/%d subset fits usable",
      elapsed, usable, length(r2))
  msg("[dom] McFadden R2 total = %.6f | sum of dominance weights = %.6f",
      max(r2, na.rm = TRUE), sum(dom, na.rm = TRUE))
  
  data.table::data.table(
    item = predictors,
    dominance_r2_m = as.numeric(dom)
  )[order(-dominance_r2_m)]
}

write_scatter_labeled <- function(data, x, y, label, out_png,
                                  xlab, ylab, title = NULL,
                                  point_size = 2.2,
                                  text_size = 3.8) {
  stopifnot(all(c(x, y, label) %in% names(data)))
  
  df <- as.data.frame(data)
  df <- df[is.finite(df[[x]]) & is.finite(df[[y]]), , drop = FALSE]
  
  xr <- range(df[[x]], na.rm = TRUE)
  yr <- range(df[[y]], na.rm = TRUE)
  xpad <- 0.06 * diff(xr)
  ypad <- 0.08 * diff(yr)
  
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey70", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey70", linewidth = 0.5) +
    ggplot2::geom_point(size = point_size) +
    ggrepel::geom_text_repel(
      ggplot2::aes_string(label = label),
      size = text_size,
      max.overlaps = Inf,
      box.padding = 0.30,
      point.padding = 0.20,
      min.segment.length = 0,
      seed = 42,
      segment.alpha = 0.5
    ) +
    ggplot2::labs(
      title = title %||% "",
      x = xlab,
      y = ylab
    ) +
    ggplot2::coord_cartesian(
      xlim = c(xr[1] - xpad, xr[2] + xpad),
      ylim = c(yr[1] - ypad, yr[2] + ypad),
      clip = "off"
    ) +
    theme_pub_local(11) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 35, 10, 10)
    )
  
  save_plot_gg_local(
    plot = p,
    filename = out_png,
    width = 8.2,
    height = 6.8,
    dpi = 300
  )
  
  invisible(p)
}

dx <- fread(DX_FILE)
psy <- fread(PSY_FILE)

if (!(ID_VAR %in% names(dx))) stop(sprintf("ID column '%s' not found in %s", ID_VAR, DX_FILE))
if (!(ID_VAR %in% names(psy))) stop(sprintf("ID column '%s' not found in %s", ID_VAR, PSY_FILE))
if (!(OUTCOME_VAR %in% names(dx))) stop(sprintf("Outcome '%s' not found in %s", OUTCOME_VAR, DX_FILE))

missing_preds <- setdiff(PREDICTORS, names(psy))
if (length(missing_preds)) stop("Missing predictors in psychometric matrix: ", paste(missing_preds, collapse = ", "))

dx[[ID_VAR]]  <- trimws(as.character(dx[[ID_VAR]]))
psy[[ID_VAR]] <- trimws(as.character(psy[[ID_VAR]]))

dx <- unique(dx[, c(ID_VAR, OUTCOME_VAR), with = FALSE], by = ID_VAR)
DT <- merge(dx, psy[, c(ID_VAR, PREDICTORS), with = FALSE], by = ID_VAR, all = FALSE)

DT[[OUTCOME_VAR]] <- coerce_binary01(DT[[OUTCOME_VAR]])
for (v in PREDICTORS) DT[[v]] <- char_to_num(DT[[v]])

keep_cc <- complete.cases(DT[, c(OUTCOME_VAR, PREDICTORS), with = FALSE])
DT <- DT[keep_cc]

y <- as.integer(DT[[OUTCOME_VAR]])
X <- as.matrix(DT[, ..PREDICTORS])

msg("[structure_std] N=%d | cases=%d | controls=%d", nrow(DT), sum(y == 1), sum(y == 0))

ridge_res <- oof_ridge_bin(X, y, K = K_FOLDS, seed = SEED, n_perm = N_PERM)
DT[, eta_oof := ridge_res$eta_oof]
DT[, p_oof   := ridge_res$p_oof]

Xz <- as.data.table(lapply(PREDICTORS, function(v) as.numeric(scale(DT[[v]]))))
setnames(Xz, PREDICTORS)

dom_dt <- NULL
dom_file <- NA_character_
dom_vars <- character(0)

if (isTRUE(DOM_ENABLED)) {
  ridge_abs <- abs(ridge_res$coef_mean[PREDICTORS])
  ridge_abs[!is.finite(ridge_abs)] <- -Inf
  
  dom_k <- min(DOM_TOP_K, length(PREDICTORS))
  dom_vars <- names(sort(ridge_abs, decreasing = TRUE))[seq_len(dom_k)]
  
  msg("[dom] selected top-%d predictors by |ridge_coef_mean|: %s",
      dom_k, paste(dom_vars, collapse = ", "))
  
  dom_dt <- dominance_mcfadden_exact(
    X = as.data.frame(Xz),
    y = y,
    predictors = dom_vars,
    ncores = DOM_CORES,
    maxit = DOM_MAXIT
  )
  
  if (!is.null(dom_dt)) {
    stopifnot(all(dom_dt$item %in% dom_vars))
    stopifnot(all(is.finite(dom_dt$dominance_r2_m)))
  }
  
  dom_file <- file.path(OUT_DIR, sprintf("dominance_top%d_r2m.csv", length(dom_vars)))
  fwrite(dom_dt, dom_file)
  msg("[dom] wrote %s", dom_file)
}

item_summary <- rbindlist(lapply(PREDICTORS, function(v) {
  xz <- Xz[[v]]
  uni <- logit_or(y, xz)
  data.table(
    item = v,
    marginal_beta = uni$beta,
    marginal_or_per_sd = uni$or,
    marginal_or_lo = uni$lo,
    marginal_or_hi = uni$hi,
    marginal_p = uni$p,
    marginal_auc = calc_auc(y, xz),
    loading_score_r = safe_cor(xz, DT$eta_oof, method = "pearson"),
    loading_score_r_spearman = safe_cor(xz, DT$eta_oof, method = "spearman"),
    ridge_coef_mean = unname(ridge_res$coef_mean[v]),
    ridge_coef_sd = unname(ridge_res$coef_sd[v]),
    ridge_sign_stability = unname(ridge_res$coef_sign_stability[v])
  )
}))

if (!is.null(ridge_res$perm_summary)) {
  item_summary <- merge(item_summary, ridge_res$perm_summary, by = "item", all.x = TRUE)
}
if (!is.null(dom_dt)) {
  item_summary <- merge(item_summary, dom_dt, by = "item", all.x = TRUE)
}

item_summary[, abs_loading := abs(loading_score_r)]
setorder(item_summary, -abs_loading)
item_summary[, abs_loading := NULL]
item_summary[, item_label := sub("^psy\\.", "", item)]

subject_scores <- DT[, c(ID_VAR, OUTCOME_VAR, "eta_oof", "p_oof"), with = FALSE]
auc_oof <- calc_auc(y, DT$p_oof)
ll_oof  <- logloss(y, DT$p_oof)
metrics_oof <- compute_binary_metrics_with_ci(y, DT$p_oof, boot_B = BOOT_B, boot_seed = BOOT_SEED)
metrics_dt <- data.table(
  metric = c('auc','auprc','auprg','logloss','brier','calibration_intercept','calibration_slope','ece'),
  point = c(metrics_oof$auc, metrics_oof$auprc, metrics_oof$auprg, metrics_oof$logloss, metrics_oof$brier, metrics_oof$cal_int, metrics_oof$cal_slope, metrics_oof$cal_ece),
  lo = c(metrics_oof$auc_lo, metrics_oof$auprc_lo, metrics_oof$auprg_lo, metrics_oof$logloss_lo, metrics_oof$brier_lo, metrics_oof$cal_int_lo, metrics_oof$cal_slope_lo, metrics_oof$cal_ece_lo),
  hi = c(metrics_oof$auc_hi, metrics_oof$auprc_hi, metrics_oof$auprg_hi, metrics_oof$logloss_hi, metrics_oof$brier_hi, metrics_oof$cal_int_hi, metrics_oof$cal_slope_hi, metrics_oof$cal_ece_hi)
)

summary_txt <- c(
  "Structure-coefficient / relative-importance analysis",
  sprintf("Diagnosis file: %s", DX_FILE),
  sprintf("Psychometric file: %s", PSY_FILE),
  sprintf("Outcome: %s", OUTCOME_VAR),
  sprintf("Predictors total: %d", length(PREDICTORS)),
  sprintf("Complete-case N = %d", nrow(DT)),
  sprintf("Cases = %d | Controls = %d", sum(y == 1), sum(y == 0)),
  sprintf("Outer CV folds = %d", K_FOLDS),
  sprintf("Permutation repetitions per predictor per outer fold = %d", N_PERM),
  sprintf("OOF AUC = %.3f [%.3f, %.3f]", metrics_oof$auc, metrics_oof$auc_lo, metrics_oof$auc_hi),
  sprintf("OOF AUPRC = %.3f [%.3f, %.3f]", metrics_oof$auprc, metrics_oof$auprc_lo, metrics_oof$auprc_hi),
  sprintf("OOF AUPRG = %.3f [%.3f, %.3f]", metrics_oof$auprg, metrics_oof$auprg_lo, metrics_oof$auprg_hi),
  sprintf("OOF LogLoss = %.3f [%.3f, %.3f]", metrics_oof$logloss, metrics_oof$logloss_lo, metrics_oof$logloss_hi),
  sprintf("OOF Brier = %.3f [%.3f, %.3f]", metrics_oof$brier, metrics_oof$brier_lo, metrics_oof$brier_hi),
  sprintf("Calibration intercept = %.3f [%.3f, %.3f]", metrics_oof$cal_int, metrics_oof$cal_int_lo, metrics_oof$cal_int_hi),
  sprintf("Calibration slope = %.3f [%.3f, %.3f]", metrics_oof$cal_slope, metrics_oof$cal_slope_lo, metrics_oof$cal_slope_hi),
  sprintf("ECE = %.3f [%.3f, %.3f]", metrics_oof$cal_ece, metrics_oof$cal_ece_lo, metrics_oof$cal_ece_hi),
  if (is.null(dom_dt)) {
    "Exact McFadden general dominance not computed."
  } else {
    sprintf(
      "Exact McFadden general dominance computed on top-%d predictors selected by |ridge_coef_mean|; output file: %s",
      length(dom_vars), basename(dom_file)
    )
  },
  "",
  "Interpretation:",
  "- marginal_or_per_sd: univariate item-outcome association",
  "- loading_score_r: structure coefficient (r)",
  "- ridge_coef_mean: mean ridge coefficient",
  "- perm_dAUC_mean / perm_dLogLoss_mean: held-out performance loss when that item is permuted",
  "- dominance_r2_m: exact general dominance weight for McFadden pseudo-R2 on the reduced predictor set"
)

writeLines(summary_txt, con = file.path(OUT_DIR, "summary.txt"))

fwrite(item_summary, file.path(OUT_DIR, "item_summary.csv"))
fwrite(subject_scores, file.path(OUT_DIR, "subject_scores.csv"))
fwrite(metrics_dt, file.path(OUT_DIR, "model_metrics.csv"))

write_binary_curves(
  y01 = y,
  p_hat = DT$p_oof,
  out_png = file.path(OUT_DIR, "roc_oof.png"),
  title = "OOF ridge-logit diagnostics",
  subtitle = sprintf("N = %d | cases = %d | controls = %d", nrow(DT), sum(y == 1), sum(y == 0)),
  metrics = metrics_oof
)

write_scatter_labeled(
  data = item_summary,
  x = "loading_score_r",
  y = "ridge_coef_mean",
  label = "item_label",
  out_png = file.path(OUT_DIR, "loadings_vs_weights.png"),
  xlab = "Structure coefficient (r)",
  ylab = "Mean ridge coefficient",
  title = "Structure coefficient vs mean ridge coefficient"
)

write_scatter_labeled(
  data = item_summary,
  x = "marginal_beta",
  y = "loading_score_r",
  label = "item_label",
  out_png = file.path(OUT_DIR, "marginal_vs_loading.png"),
  xlab = "Univariate log-odds coefficient (per s.d.)",
  ylab = "Structure coefficient (r)",
  title = "Univariate log-odds coefficient vs structure coefficient"
)

msg("[structure_std] wrote outputs to %s", OUT_DIR)
msg("[structure_std] OOF AUC = %.3f", auc_oof)
