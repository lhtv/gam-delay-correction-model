## ---- Setup ----
library(data.table)
library(mgcv)
library(ggplot2)
library(mvtnorm)
library(readr)
library(gridExtra)
library(dplyr)
library(knitr)
library(kableExtra)

set.seed(123)

setwd("D:/Dissertation")
plotDir <- file.path(getwd(), "Plots")
if (!dir.exists(plotDir)) dir.create(plotDir, recursive = TRUE)

## ---- Load data (only once) ----
if (!exists("dengue.data") || !exists("dengue.data.full")) {
  load("dengue2013April15thEpiWeek.RData")
  load("dengue2013April15thEpiWeekFull.RData")
}
dengue.data      <- data.table(dengue.data)
dengue.data.full <- data.table(dengue.data.full)

## ---- Basic cleaning and column setup ----
if ("Obs" %in% names(dengue.data.full)) dengue.data.full[, Obs := NULL]
index <- -c(1:5)
names(dengue.data.full)[index] <- paste0("delay", names(dengue.data.full)[index])
setnames(dengue.data.full, "Total", "y")

## ---- Choose delays and derive x ----
Ndelay <- 10
delay_index <- grep("^delay", names(dengue.data.full))
stopifnot(length(delay_index) >= (Ndelay + 1))

keep_cols <- c(1:5, delay_index[1:(Ndelay + 1)])
dat <- dengue.data.full[, ..keep_cols]
delay_indexHere <- delay_index[1:(Ndelay + 1)]
dat[, x := rowSums(.SD), .SDcols = names(dat)[delay_indexHere]]

## ---- Evaluation weeks ----
k_seq <- 160:170

## ---- Storage objects ----
width_tab   <- list()  # PI width rows across k and models
acc_tab     <- list()  # accuracy rows across k and models
metrics_tab <- list()  # per-k metric tables for both models

## ---- Helper: mask function ----
mask_for_k <- function(DT, k, delay_idx_here, Ndelay) {
  out <- copy(DT)
  delay_cols <- names(out)[delay_idx_here]
  
  # Nowcasting mask: hide "next" delay column for the last Ndelay weeks
  for (i in 1:Ndelay) {
    ind <- delay_idx_here[i] + 1L
    if (ind <= ncol(out)) {
      colname <- names(out)[ind]
      r1 <- max(1L, k - i + 1L)
      r2 <- k
      if (r1 <= r2) out[r1:r2, (colname) := NA_real_]
    }
  }
  
  # Forecasting mask: hide all delay columns after week k
  if (k + 1L <= nrow(out)) {
    out[(k + 1L):.N, (delay_cols) := NA_real_]
  }
  
  # Recompute x after masking
  out[, x := rowSums(.SD), .SDcols = delay_cols]
  out
}

## ---- Metrics helpers ----
rmse     <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae      <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
cov95    <- function(obs, lo, hi) mean(obs >= lo & obs <= hi, na.rm = TRUE)
pi_width <- function(lo, hi) mean(hi - lo, na.rm = TRUE)

## ---- Simulation settings ----
n.sims <- 1000


## Main loop over k

for (k in k_seq) {
  
  ## ----- Mask data for this k -----
  datS <- mask_for_k(dat, k, delay_indexHere, Ndelay)
  
  
  ## Model XZ: joint gfam for delays and x
  
  longDat_XZ <- melt(
    datS,
    id.vars = c("DATE", "EpiWeek", "Time", "x"),
    measure.vars = c(grep("^delay", colnames(datS), value = TRUE), "x"),
    variable.name = "variable", value.name = "value"
  )
  longDat_XZ[, Nvariable := as.integer(factor(variable, levels = unique(variable)))]
  longDat_XZ[, omega := factor(Nvariable)]
  longDat_XZ[, y_resp := value]
  
  longDf_XZ <- as.data.frame(longDat_XZ)
  longDf_XZ$y <- cbind(longDf_XZ$y_resp, longDf_XZ$Nvariable)
  
  n_comp_XZ <- length(unique(longDf_XZ$Nvariable))
  
  XZmodel <- gam(
    y ~ omega +
      s(Time, k = 20, bs = "tp", m = 1) +
      s(Time, omega, bs = "sz", xt = list(bs = "tp"), k = 5, m = 1) +
      s(EpiWeek, k = 15, bs = "cc"),
    family = gfam(rep(list(nb), n_comp_XZ)),
    method = "REML",
    data = longDf_XZ,
    knots = list(EpiWeek = c(0, 52))
  )
  
  ## Predictive simulation for XZ (all components)
  predDat_XZ <- data.table(longDat_XZ)
  predDat_XZ[, y := Nvariable]  # gfam expects y = component index for prediction
  
  X_lp_XZ <- predict(XZmodel, newdata = predDat_XZ, type = "lpmatrix")
  b_XZ    <- mvtnorm::rmvnorm(n.sims, mean = coef(XZmodel), sigma = XZmodel$Vp)
  LP_XZ   <- tcrossprod(X_lp_XZ, b_XZ)
  Mean_XZ <- exp(LP_XZ)
  
  thetas_XZ <- exp(XZmodel$family$getTheta())
  if (length(thetas_XZ) < n_comp_XZ) thetas_XZ <- rep(thetas_XZ[1], n_comp_XZ)
  theta_by_row_XZ <- thetas_XZ[predDat_XZ$Nvariable]
  
  Y_XZ <- matrix(NA_real_, nrow = nrow(Mean_XZ), ncol = ncol(Mean_XZ))
  for (j in 1:ncol(Mean_XZ)) {
    Y_XZ[, j] <- rnbinom(nrow(Mean_XZ), mu = Mean_XZ[, j], size = theta_by_row_XZ)
  }
  
  predDat_XZ[, preds := apply(Y_XZ, 1, stats::median)]
  predDat_XZ[, lwr   := apply(Y_XZ, 1, stats::quantile, probs = 0.025)]
  predDat_XZ[, upr   := apply(Y_XZ, 1, stats::quantile, probs = 0.975)]
  
  x_level_XZ <- unique(longDat_XZ[variable == "x", Nvariable])[1]
  width_XZ <- predDat_XZ[Nvariable == x_level_XZ,
                         .(Time, width = pmax(upr - lwr, 0))]
  width_XZ[, `:=`(k = k, model = "XZ")]
  
  ## For accuracy metrics, attach DATE and obs for x component
  predXZ_x <- predDat_XZ[Nvariable == x_level_XZ, .(DATE, Time, preds, lwr, upr)]
  predXZ_x <- merge(predXZ_x, dat[, .(DATE, y)], by = "DATE", all.x = TRUE)
  
  
  ## Model Z: GAM on delays, then sum to x
  
  longDat_Z <- melt(
    datS,
    id.vars = c("DATE", "EpiWeek", "Time"),
    measure.vars = grep("^delay", colnames(datS), value = TRUE),
    variable.name = "variable", value.name = "value"
  )
  longDat_Z[, delay := as.integer(gsub("delay", "", variable))]
  
  Zmodel <- gam(
    value ~ s(Time, k = 20, bs = "tp") +
      s(delay, k = 8, bs = "tp") +
      ti(Time, delay, bs = "tp", k = 5) +
      s(EpiWeek, k = 15, bs = "cc"),
    family = nb,
    method = "REML",
    data = longDat_Z,
    knots = list(EpiWeek = c(0, 52))
  )
  
  ## Predictive simulation for Z, then aggregate to x by DATE
  X_lp_Z <- predict(Zmodel, newdata = longDat_Z, type = "lpmatrix")
  b_Z    <- mvtnorm::rmvnorm(n.sims, mean = coef(Zmodel), sigma = Zmodel$Vp)
  LP_Z   <- tcrossprod(X_lp_Z, b_Z)
  Mean_Z <- exp(LP_Z)
  
  thetas_Z <- exp(Zmodel$family$getTheta())
  Y_Z <- apply(Mean_Z, 2, function(mu) rnbinom(length(mu), mu = mu, size = thetas_Z))
  
  longDatPred_Z  <- cbind(longDat_Z, Y_Z)
  shortDatPred_Z <- longDatPred_Z[, lapply(.SD, sum), by = DATE, .SDcols = paste0("V", 1:n.sims)]
  
  shortDatPred_Z[, preds := apply(.SD, 1, median), .SDcols = paste0("V", 1:n.sims)]
  shortDatPred_Z[, lwr   := apply(.SD, 1, quantile, probs = 0.025), .SDcols = paste0("V", 1:n.sims)]
  shortDatPred_Z[, upr   := apply(.SD, 1, quantile, probs = 0.975), .SDcols = paste0("V", 1:n.sims)]
  
  pred_Z <- merge(dat[, .(DATE, Time, y)],
                  shortDatPred_Z[, .(DATE, preds, lwr, upr)],
                  by = "DATE", all.x = TRUE)
  
  width_Z <- pred_Z[, .(Time, width = pmax(upr - lwr, 0))]
  width_Z[, `:=`(k = k, model = "Z")]
  
  
  ## k-check
  k.check(Zmodel)
  k.check(XZmodel)
  
  ## QQ plots
  qq.gam(Zmodel)
  qq.gam(XZmodel)
  
  ## Collect PI widths and accuracy rows
  
  width_tab[[length(width_tab) + 1L]] <- rbindlist(list(width_XZ, width_Z), use.names = TRUE)
  
  acc_XZ_k <- predXZ_x[, .(Time, k, model = "XZ", obs = y, preds, lwr, upr)]
  acc_Z_k  <- pred_Z[,   .(Time, k, model = "Z",  obs = y, preds, lwr, upr)]
  acc_tab[[length(acc_tab) + 1L]] <- rbindlist(list(acc_XZ_k, acc_Z_k), use.names = TRUE)
  
  
  ## Per-k metrics and store
  
  metrics_XZ_k <- data.frame(
    k = k,
    Model = "Model XZ",
    AIC = AIC(XZmodel),
    `Deviance explained` = summary(XZmodel)$dev.expl,
    RMSE = rmse(predXZ_x$y, predXZ_x$preds),
    MAE  = mae (predXZ_x$y, predXZ_x$preds),
    `95% IC` = cov95(predXZ_x$y, predXZ_x$lwr, predXZ_x$upr),
    `PI Width` = pi_width(predXZ_x$lwr, predXZ_x$upr),
    check.names = FALSE
  )
  
  metrics_Z_k <- data.frame(
    k = k,
    Model = "Model Z",
    AIC = AIC(Zmodel),
    `Deviance explained` = summary(Zmodel)$dev.expl,
    RMSE = rmse(pred_Z$y, pred_Z$preds),
    MAE  = mae (pred_Z$y, pred_Z$preds),
    `95% IC` = cov95(pred_Z$y, pred_Z$lwr, pred_Z$upr),
    `PI Width` = pi_width(pred_Z$lwr, pred_Z$upr),
    check.names = FALSE
  )
  
  metrics_tab[[length(metrics_tab) + 1L]] <- rbind(metrics_Z_k, metrics_XZ_k)
}


## Combine results across k

width_all <- rbindlist(width_tab, use.names = TRUE)
acc_all   <- rbindlist(acc_tab, use.names = TRUE)

## Horizons per row, based on each row's k
width_all[, horizon := ifelse(Time <= k, "Nowcast", "Forecast")]
width_all[, horizon := factor(horizon, levels = c("Nowcast", "Forecast"))]

acc_all[, horizon := ifelse(Time <= k, "Nowcast", "Forecast")]
acc_all[, abs_err := abs(preds - obs)]
acc_all[, sq_err  := (preds - obs)^2]


## Plots

# PI width by k and model
# p1 <- ggplot(width_all, aes(x = factor(k), y = width, fill = model)) +
#   theme_bw() +
#   geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8)) +
#   labs(
#     x = "Evaluation week k",
#     y = "Predictive interval width",
#     title = "PI width comparison by model across evaluation weeks",
#     fill = "Model"
#   )

# Overall PI width by model
p2 <- ggplot(width_all, aes(x = model, y = width, fill = model)) +
  geom_boxplot(width = 0.5) +
  theme_bw() +
  labs(
    x = "Model",
    y = "Predictive interval width",
    title = "Overall PI width comparison between models"
  ) +
  theme(legend.position = "none")

# PI width by horizon and model
p3 <- ggplot(width_all, aes(x = horizon, y = width, fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  theme_bw() +
  labs(
    x = "Horizon",
    y = "Predictive interval width",
    title = "PI width comparison (160 ≤ k ≤ 170) by model and horizon"
  )

# Absolute error by horizon and model
p4 <- ggplot(acc_all, aes(x = horizon, y = abs_err, fill = model)) +
  geom_boxplot(
    aes(group = interaction(horizon, model)),
    position = position_dodge2(width = 0.75, preserve = "single"),
    outlier.size = 1.2
  ) +
  stat_summary(
    aes(color = model, group = model),
    fun = median, geom = "point",
    position = position_dodge2(width = 0.75, preserve = "single"),
    size = 3, shape = 18
  ) +
  theme_bw() +
  labs(
    x = "Horizon",
    y = "Absolute error",
    title = "Prediction accuracy (160 ≤ k ≤ 170): absolute error by model and horizon"
  )

# Print and optionally save
print(p1); print(p2); print(p3); print(p4)
ggsave(file.path(plotDir, "pi_width_by_k_model.png"), p1, width = 9, height = 5, dpi = 300)
ggsave(file.path(plotDir, "pi_width_overall.png"),     p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(plotDir, "pi_width_by_horizon.png"),  p3, width = 7, height = 5, dpi = 300)
ggsave(file.path(plotDir, "abs_error_by_horizon.png"), p4, width = 7, height = 5, dpi = 300)


## Metrics tables: per-k and averaged

metrics_all <- bind_rows(metrics_tab)

# # Per-k table
# metrics_all %>%
#   arrange(k, Model) %>%
#   kable(format = "latex", booktabs = TRUE, digits = 3,
#         caption = "Model performance by evaluation week k") %>%
#   kable_styling(latex_options = c("hold_position")) %>%
#   cat(file = "performance_comparison_by_k.tex")

# Averaged across k
metrics_mean <- metrics_all %>%
  group_by(Model) %>%
  summarise(
    AIC = mean(AIC, na.rm = TRUE),
    `Deviance explained` = mean(`Deviance explained`, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    MAE  = mean(MAE,  na.rm = TRUE),
    `95% IC` = mean(`95% IC`, na.rm = TRUE),
    `PI Width` = mean(`PI Width`, na.rm = TRUE),
    .groups = "drop"
  )

metrics_mean %>%
  kable(format = "latex", booktabs = TRUE, digits = 3,
        caption = "Average model performance across k = 160–170") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  cat(file = "performance_comparison_average.tex")
