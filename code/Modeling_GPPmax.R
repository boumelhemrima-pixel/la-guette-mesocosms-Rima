# =============================
#      Load Required Libraries
# =============================
library(hydroGOF)
library(dplyr)
library(minpack.lm)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(scales)
library(tidyr)
library(tibble)
library(readxl)

# =============================
#        Load and Prepare Data
# =============================
cbd <- read.csv("C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by campaign/by replicate/mesocosm_data1.csv", sep = ",")
cbd <- dplyr::rename(cbd, Date = TIMESTAMP.Reco)
cbd$Date <- as.POSIXct(cbd$Date, format = "%m/%d/%Y %H:%M")
GPPmax_cal <- cbd

Tmin <- 0; Tmax <- 45

setwd("C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by campaign/by replicate/GPP_by treatment/for git")

veg_types <- c("B", "BG", "BGS", "GS", "G")
num_iter <- 1
results_all <- data.frame()
fitted_models <- list()
params_list <- list()
GPPmax_predictions <- data.frame()

# Desired plot orders
subset_levels_order <- c("B","BG","BGS","GS","G")
group_levels_order  <- c("B","BG","BGS","GS","G")

# === master switch: keep code, but toggle VI behavior ===
USE_VI <- FALSE   # set TRUE to re-enable VI everywhere

# =============================
#        Fit Models and Compute Residuals
# =============================

bounded_temp <- function(Ta, Topt, Tmin = 0, Tmax = 45) {
  num <- (Ta - Tmin) * (Ta - Tmax)
  den <- num - (Ta - Topt)^2
  return(ifelse(den == 0, NA, num / den))
}

fit_models_for_group <- function(veg_name, GPPmax_cal, num_iter = 1) {
  # Subset all replicates for the given vegetation group
  veg_data <- subset(GPPmax_cal, Plot == veg_name)
  
  # Apply filters once per group (merged replicates)
  veg_data <- veg_data %>%
    filter(!is.na(GPPmax), !is.na(LAI), (!USE_VI | !is.na(VI)), !is.na(Ta), GPPmax < 0) %>%
    mutate(across(c(GPPmax, VI, LAI, Ta), as.numeric))
  
  for (iteration in 1:num_iter) {
    for (model_type in (if (USE_VI) c("Model_VI", "Model_LAI") else c("Model_LAI"))) {
      tryCatch({
        # Define model equation
        formula_model <- switch(model_type,
                                "Model_VI"  = GPPmax ~ ((a * VI  + b) * (((Ta - Tmin) * (Ta - Tmax)) / (((Ta - Tmin) * (Ta - Tmax)) - (Ta - Topt)^2))),
                                "Model_LAI" = GPPmax ~ ((a * LAI + b) * (((Ta - Tmin) * (Ta - Tmax)) / (((Ta - Tmin) * (Ta - Tmax)) - (Ta - Topt)^2)))
        )
        
        # Fit per vegetation (not replicate)
        fit <- nlsLM(
          formula_model,
          data   = veg_data,
          start  = list(a = 1, b = 1, Topt = 10),
          lower  = c(a = -Inf, b = -Inf, Topt = 0),
          upper  = c(a =  Inf, b =  Inf, Topt = Inf),
          control = nls.lm.control(maxiter = 100000, ftol = 1e-10)
        )
        
        preds <- predict(fit)
        adj_r2_val <- summary(lm(veg_data$GPPmax ~ preds))$adj.r.squared
        rmse_val   <- rmse(veg_data$GPPmax, preds)
        aic_val    <- AIC(fit)
        bic_val    <- BIC(fit)
        
        model_key <- paste0(veg_name, "_iter", iteration, "_", model_type)
        fitted_models[[model_key]] <<- fit
        params_list[[model_key]]  <<- coef(fit)
        
        s <- summary(fit)
        coefs   <- coef(s)
        errors  <- as.numeric(coefs[, "Std. Error"])
        params  <- as.numeric(coefs[, "Estimate"])
        t_stat  <- params / errors
        df_fit  <- df.residual(fit)
        p_vals  <- 2 * pt(-abs(t_stat), df = df_fit)
        
        row <- list(
          SubsetType = veg_name, Iteration = iteration, Model = model_type,
          AIC = aic_val, BIC = bic_val, RMSE = rmse_val, R2 = adj_r2_val,
          Veg = veg_name,
          a = params[1], b = params[2], Topt = params[3],
          SD_a = errors[1], SD_b = errors[2], SD_Topt = errors[3],
          pval_a = p_vals[1], pval_b = p_vals[2], pval_Topt = p_vals[3]
        )
        
        results_all <<- bind_rows(results_all, as.data.frame(row))
        
        GPPmax_predictions <<- bind_rows(GPPmax_predictions, data.frame(
          Veg = veg_name, Model = model_type,
          Date = veg_data$Date,
          GPPmax_obs = veg_data$GPPmax,
          GPPmax_mod = preds,
          Residual = veg_data$GPPmax - preds
        ))
        
      }, error = function(e) {
        cat("❌ Error in", veg_name, model_type, ":", e$message, "\n")
      })
    }
  }
}

# --- Run for each vegetation group
for (veg in veg_types) {
  fit_models_for_group(veg, GPPmax_cal, num_iter)
}

# =============================
# Replace Date column in GPPmax_predictions with original timestamp
# =============================
cbd_with_ts <- cbd %>%
  rename(`TIMESTAMP.Reco` = Date) %>%
  mutate(Plot_ID = Plot)

GPPmax_predictions <- GPPmax_predictions %>%
  mutate(Plot_ID = Veg) %>%
  left_join(
    cbd_with_ts %>%
      select(Plot_ID, `TIMESTAMP.Reco`, GPPmax) %>%
      rename(GPPmax_obs_check = GPPmax),
    by = c("Plot_ID", "GPPmax_obs" = "GPPmax_obs_check")
  ) %>%
  select(-Date) %>%
  rename(Date = `TIMESTAMP.Reco`) %>%
  select(-Plot_ID)

GPPmax_predictions$Date <- format(GPPmax_predictions$Date, "%m/%d/%Y %H:%M")

# =============================
# Save results to Excel
# =============================
write.xlsx(results_all,        file = "GPPmax_Model_Results_All.xlsx",           rowNames = FALSE)
write.xlsx(GPPmax_predictions, file = "GPPmax_Model_Residuals_with_Scaled.xlsx", rowNames = FALSE)

# =============================================================
#   EXPAND GROUP PARAMETERS → REPLICATES + COMPUTE R² / RMSE / p-value
# =============================================================

library(openxlsx)
library(dplyr)
library(hydroGOF)

# --- load model results
res_path <- "GPPmax_Model_Results_All.xlsx"
if (!file.exists(res_path)) stop("❌  GPPmax_Model_Results_All.xlsx not found.")
res_df <- read.xlsx(res_path)

# --- mapping Veg → replicates
expand_map <- list(
  B   = c("B1","B2","B3"),
  BG  = c("BG1","BG2","BG3"),
  BGS = c("BGS1","BGS2","BGS3"),
  GS  = c("GS1","GS2","GS3"),
  G   = c("G1","G2","G3")
)

# --- helper: bounded-temperature function (same as in your model)
bounded_temp <- function(Ta, Topt, Tmin = 0, Tmax = 45) {
  num <- (Ta - Tmin) * (Ta - Tmax)
  den <- num - (Ta - Topt)^2
  val <- num / den
  val[!is.finite(val)] <- NA
  val
}

# --- compute R², RMSE, p-value for one replicate
calc_metrics <- function(rep_data, a, b, Topt) {
  if (nrow(rep_data) < 5) return(c(R2 = NA, RMSE = NA, p_value = NA))
  preds <- (a * rep_data$LAI + b) * bounded_temp(rep_data$Ta, Topt)
  preds[!is.finite(preds)] <- NA
  keep <- is.finite(preds) & is.finite(rep_data$GPPmax)
  if (sum(keep) < 3) return(c(R2 = NA, RMSE = NA, p_value = NA))
  lm_fit <- lm(rep_data$GPPmax[keep] ~ preds[keep])
  R2  <- summary(lm_fit)$adj.r.squared
  RMSE <- hydroGOF::rmse(rep_data$GPPmax[keep], preds[keep])
  pval <- summary(lm_fit)$coefficients[2,4]
  c(R2 = R2, RMSE = RMSE, p_value = pval)
}

# --- main loop: duplicate and compute
expanded_list <- lapply(seq_len(nrow(res_df)), function(i) {
  row <- res_df[i,]
  veg <- as.character(row$Veg)
  if (!veg %in% names(expand_map)) return(NULL)
  
  reps <- expand_map[[veg]]
  # extract parameters
  a <- row$a; b <- row$b; Topt <- row$Topt
  
  do.call(rbind, lapply(reps, function(rid) {
    rep_id <- sub(veg, "", rid)
    rep_data <- subset(GPPmax_cal, Plot == veg & Replicates == rep_id)
    m <- calc_metrics(rep_data, a, b, Topt)
    out <- row
    out$SubsetType <- rid
    out$Replicate  <- rep_id
    out$R2_repl    <- m["R2"]
    out$RMSE_repl  <- m["RMSE"]
    out$pval_repl  <- m["p_value"]
    out
  }))
})

expanded_df <- bind_rows(expanded_list)

# reorder columns
expanded_df <- expanded_df %>%
  select(SubsetType, Replicate, Veg, Model, everything(),
         R2_repl, RMSE_repl, pval_repl)

# save new version
write.xlsx(expanded_df, "GPPmax_Model_Results_All_withReplicates.xlsx", rowNames = FALSE)
cat("✅  Saved: GPPmax_Model_Results_All_withReplicates.xlsx  (duplicated params + R²/RMSE/p)\n")

# =============================
# Merge predictions back to cbd (wide per model)
# =============================
residuals_long <- GPPmax_predictions %>%
  mutate(Date = as.POSIXct(Date, format = "%m/%d/%Y %H:%M"))

residuals_summarized <- residuals_long %>%
  group_by(Date, Veg, Model) %>%
  summarise(
    GPPmax_obs = mean(GPPmax_obs, na.rm = TRUE),
    GPPmax_mod = mean(GPPmax_mod, na.rm = TRUE),
    Residual   = mean(Residual, na.rm = TRUE),
    .groups = "drop"
  )

residuals_wide <- residuals_summarized %>%
  pivot_wider(
    id_cols = c(Date, Veg),
    names_from = Model,
    values_from = c(GPPmax_obs, GPPmax_mod, Residual),
    names_glue = "{.value}_{Model}"
  ) %>%
  rename(Plot = Veg)

cbd_with_all_residuals <- cbd %>%
  left_join(residuals_wide, by = c("Date", "Plot"))

write.csv(cbd_with_all_residuals, "mesocosm_data1_with_residuals.csv", row.names = FALSE)
cat("✅ Residuals and modeled GPPmax values saved into mesocosm_data1_with_residuals.csv\n")

# =============================
#        Performance Plotting (per-fit scatter)
# =============================
# =============================
#   Performance Plotting (per-fit scatter) — replicate level
# =============================

# --- Expand group-level parameters to all replicates
expand_map <- list(
  B   = c("B1","B2","B3"),
  BG  = c("BG1","BG2","BG3"),
  BGS = c("BGS1","BGS2","BGS3"),
  GS  = c("GS1","GS2","GS3"),
  G   = c("G1","G2","G3")
)

# Expand results_all to replicate level
results_repl <- do.call(rbind, lapply(seq_len(nrow(results_all)), function(i) {
  veg <- as.character(results_all$Veg[i])
  if (veg %in% names(expand_map)) {
    reps <- expand_map[[veg]]
    out <- results_all[rep(i, length(reps)), ]
    out$SubsetType <- reps
    out$Replicate <- sub(veg, "", reps)
    return(out)
  } else {
    return(results_all[i, , drop = FALSE])
  }
}))

# --- Recompute Adjusted R² and RMSE per replicate using duplicated parameters
compute_metrics <- function(df, a, b, Topt) {
  # Only valid if required columns exist
  if (all(c("LAI","Ta","GPPmax") %in% names(df))) {
    preds <- (a * df$LAI + b) *
      (((df$Ta - Tmin) * (df$Ta - Tmax)) /
         (((df$Ta - Tmin) * (df$Ta - Tmax)) - (df$Ta - Topt)^2))
    preds[!is.finite(preds)] <- NA_real_
    adjR2 <- tryCatch(summary(lm(df$GPPmax ~ preds))$adj.r.squared, error = function(e) NA)
    RMSE  <- tryCatch(hydroGOF::rmse(df$GPPmax, preds), error = function(e) NA)
    return(c(R2 = adjR2, RMSE = RMSE))
  } else {
    return(c(R2 = NA, RMSE = NA))
  }
}

# --- Apply per replicate
metrics_list <- lapply(seq_len(nrow(results_repl)), function(i) {
  row <- results_repl[i, ]
  veg_name <- row$Veg
  subset_name <- row$SubsetType
  rep_data <- subset(GPPmax_cal, Plot == veg_name & Replicates %in% sub(veg_name, "", subset_name))
  met <- compute_metrics(rep_data, row$a, row$b, row$Topt)
  data.frame(
    SubsetType = subset_name,
    Veg = veg_name,
    Model = row$Model,
    R2 = met["R2"],
    RMSE = met["RMSE"]
  )
})

metrics_df <- do.call(rbind, metrics_list)

# --- Plot replicate-level scatter (RMSE vs R²)
shape_vals <- c("B" = 15, "BG" = 16, "BGS" = 17, "GS" = 18, "G" = 8)
p <- ggplot(metrics_df, aes(x = R2, y = RMSE, label = SubsetType, shape = Veg)) +
  geom_point(color = "black", size = 4) +
  geom_text_repel(size = 5, max.overlaps = 100, box.padding = 0.3, segment.color = "gray60") +
  facet_wrap(~ Model) +
  scale_shape_manual(values = shape_vals) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  theme_bw(base_size = 20) +
  labs(title = "Calibration RMSE vs Adjusted R² per Replicate",
       x = expression(Adjusted~R^2), y = "RMSE", shape = "Vegetation") +
  theme(strip.text = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "gray85"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18))
ggsave("Calibration_RMSE_vs_R2_by_Model_REPLICATE.png", plot = p, width = 12, height = 8, dpi = 300)
cat("✅ Saved: Calibration_RMSE_vs_R2_by_Model_REPLICATE.png (replicate-level)\n")

# =============================
#        Parameter Bar Plots (per vegetation) — with asterisks
# =============================
param_long <- results_all %>%
  pivot_longer(cols = c(a, b, Topt), names_to = "Parameter", values_to = "Estimate") %>%
  left_join(results_all %>%
              pivot_longer(cols = c(SD_a, SD_b, SD_Topt), names_to = "SD_Name", values_to = "SE") %>%
              mutate(Parameter = gsub("SD_", "", SD_Name)),
            by = c("Veg", "Model", "Parameter")) %>%
  left_join(results_all %>%
              pivot_longer(cols = c(pval_a, pval_b, pval_Topt), names_to = "Pval_Name", values_to = "Pval") %>%
              mutate(Parameter = gsub("pval_", "", Pval_Name)),
            by = c("Veg", "Model", "Parameter")) %>%
  mutate(
    Signif = ifelse(Pval < 0.05, "*", ""),
    Veg = factor(Veg, levels = subset_levels_order, ordered = TRUE),
    Model = factor(Model, levels = c("Model_VI","Model_LAI"))
  )

param_a    <- dplyr::filter(param_long, Parameter == "a")
param_b    <- dplyr::filter(param_long, Parameter == "b")
param_Topt <- dplyr::filter(param_long, Parameter == "Topt")

pdf("Parameter_a_b_Topt.pdf", width = 14, height = 8)

p_a <- ggplot(param_a, aes(x = Veg, y = Estimate)) +
  geom_col(color = "black", fill = NA, width = 0.7) +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2) +
  geom_text(aes(label = Signif, y = Estimate + SE + 0.1), vjust = 0, size = 5) +
  facet_wrap(~Model, ncol = 2, scales = "free_y") +
  labs(title = "Parameter a", y = "Estimate ± SE", x = "Vegetation") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(size = 12, face = "bold"))
print(p_a)

p_b <- ggplot(param_b, aes(x = Veg, y = Estimate)) +
  geom_col(color = "black", fill = NA, width = 0.7) +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2) +
  geom_text(aes(label = Signif, y = Estimate + SE + 0.1), vjust = 0, size = 5) +
  facet_wrap(~Model, ncol = 2, scales = "free_y") +
  labs(title = "Parameter b", y = "Estimate ± SE", x = "Vegetation") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(size = 12, face = "bold"))
print(p_b)

p_Topt <- ggplot(param_Topt, aes(x = Veg, y = Estimate)) +
  geom_col(color = "black", fill = NA, width = 0.7) +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2) +
  geom_text(aes(label = Signif, y = Estimate + SE + 0.1), vjust = 0, size = 5) +
  facet_wrap(~Model, ncol = 2, scales = "free_y") +
  labs(title = "Parameter Topt", y = "Estimate ± SE", x = "Vegetation") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(size = 12, face = "bold"))
print(p_Topt)

dev.off()

# =============================
# Calculate GPPmax predictions for entire dataset using best params
# =============================

Tmin_val <- Tmin
Tmax_val <- Tmax

bounded_temp_calc <- function(Ta, Topt) {
  num <- (Ta - Tmin_val) * (Ta - Tmax_val)
  den <- num - (Ta - Topt)^2
  out <- num / den
  out[den == 0] <- NA_real_
  out[!is.finite(out)] <- NA_real_
  out
}

# Extract best-fit parameters per vegetation and model type
params_best <- results_all %>%
  group_by(Veg, Model) %>%
  slice_min(AIC, with_ties = FALSE) %>%
  ungroup() %>%
  select(Veg, Model, a, b, Topt)

# Merge parameters back to main dataset
cbd_with_params <- cbd %>%
  left_join(
    params_best %>%
      filter(Model == "Model_VI") %>%
      rename(a_old = a, b_old = b, Topt_old = Topt),
    by = c("Plot" = "Veg")
  ) %>%
  left_join(
    params_best %>%
      filter(Model == "Model_LAI") %>%
      rename(a_new = a, b_new = b, Topt_new = Topt),
    by = c("Plot" = "Veg")
  ) %>%
  mutate(
    Ta  = as.numeric(Ta),
    VI  = as.numeric(VI),
    LAI = as.numeric(LAI),
    a_old = as.numeric(a_old), b_old = as.numeric(b_old), Topt_old = as.numeric(Topt_old),
    a_new = as.numeric(a_new), b_new = as.numeric(b_new), Topt_new = as.numeric(Topt_new)
  ) %>%
  mutate(
    GPPmax_VI_pred = ifelse(
      USE_VI &
        !is.na(a_old) & !is.na(b_old) & !is.na(Topt_old) & !is.na(VI) & !is.na(Ta),
      (a_old * VI + b_old) * bounded_temp_calc(Ta, Topt_old),
      NA_real_
    ),
    GPPmax_LAI_pred = ifelse(
      !is.na(a_new) & !is.na(b_new) & !is.na(Topt_new) & !is.na(LAI) & !is.na(Ta),
      (a_new * LAI + b_new) * bounded_temp_calc(Ta, Topt_new),
      NA_real_
    )
  )

write.csv(cbd_with_params, "mesocosm_data1_with_GPPmax_predictions.csv", row.names = FALSE)
cat(paste0("✅ GPPmax predictions (", if (USE_VI) "VI & LAI" else "LAI only", ") computed and saved.\n"))

# =============================
#   GPPmax — Averaged Parameter Plots (free y, star if ANY vegetation significant)
# =============================

gpp_long_params <- results_all %>%
  pivot_longer(cols = c(a, b, Topt), names_to = "Parameter", values_to = "Estimate") %>%
  pivot_longer(cols = c(SD_a, SD_b, SD_Topt), names_to = "SE_name", values_to = "SE") %>%
  pivot_longer(cols = c(pval_a, pval_b, pval_Topt), names_to = "Pval_name", values_to = "Pval") %>%
  mutate(
    SE_param   = sub("^SD_", "", SE_name),
    Pval_param = sub("^pval_", "", Pval_name)
  ) %>%
  filter(Parameter == SE_param, Parameter == Pval_param) %>%
  select(Veg, Model, Parameter, Estimate, SE, Pval)

gpp_summary_params <- gpp_long_params %>%
  group_by(Veg, Model, Parameter) %>%
  summarise(
    Mean_Estimate = mean(Estimate, na.rm = TRUE),
    SE_Estimate   = sd(Estimate, na.rm = TRUE) / sqrt(sum(!is.na(Estimate))),
    Min_Pval      = suppressWarnings(min(Pval, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Signif = case_when(
      is.finite(Min_Pval) & Min_Pval < 0.001 ~ "***",
      is.finite(Min_Pval) & Min_Pval < 0.01  ~ "**",
      is.finite(Min_Pval) & Min_Pval < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Model = factor(Model, levels = c("Model_VI", "Model_LAI")),
    Veg   = factor(Veg, levels = group_levels_order, ordered = TRUE)
  )

plot_param_summary_gpp <- function(param_letter) {
  df <- dplyr::filter(gpp_summary_params, Parameter == param_letter) %>% droplevels()
  panel_max <- df %>%
    group_by(Model) %>%
    summarise(ymax = max(Mean_Estimate + SE_Estimate, na.rm = TRUE), .groups = "drop")
  df <- df %>%
    left_join(panel_max, by = "Model") %>%
    mutate(
      .yhead = case_when(
        is.finite(ymax) & ymax >= 0 ~ ymax * 1.05,
        is.finite(ymax) & ymax <  0 ~ ymax * 0.95,
        TRUE ~ NA_real_
      ),
      .ylab  = case_when(
        is.finite(ymax) & ymax != 0 ~ (Mean_Estimate + SE_Estimate) + 0.02 * abs(ymax),
        TRUE ~ (Mean_Estimate + SE_Estimate) + 0.05
      )
    )
  
  ggplot(df, aes(x = Veg, y = Mean_Estimate)) +
    geom_col(fill = "white", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = Mean_Estimate - SE_Estimate,
                      ymax = Mean_Estimate + SE_Estimate), width = 0.2) +
    geom_blank(aes(y = .yhead)) +
    geom_text(aes(label = Signif, y = .ylab), size = 5, vjust = 0, na.rm = TRUE) +
    facet_wrap(~ Model, ncol = 2, scales = "free_y", drop = TRUE) +
    labs(title = paste("Parameter", param_letter),
         x = "Vegetation Group", y = "Estimate ± SE") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold", size = 16),
          strip.text  = element_text(face = "bold", size = 12))
}

pdf("GPPmax_Averaged_Parameter_Histograms.pdf", width = 12, height = 6)
for (p in c("a", "b", "Topt")) print(plot_param_summary_gpp(p))
dev.off()
cat("✅ Saved: GPPmax_Averaged_Parameter_Histograms.pdf (*** / ** / * based on min p-value per vegetation group)\n")

# ============================================================
#              PERFORMANCE TESTS & OVERALL METRICS
#              (Calibration-style, across vegetation groups)
# ============================================================

adj_r2_from_vectors <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) < 2) return(NA_real_)
  summary(lm(y[keep] ~ yhat[keep]))$adj.r.squared
}
matrix_to_long <- function(mat) {
  if (is.null(mat)) return(tibble(Model1=character(), Model2=character(), p_adj=numeric()))
  as_tibble(mat, rownames = "Model1") |>
    pivot_longer(-Model1, names_to = "Model2", values_to = "p_adj") |>
    filter(!is.na(p_adj))
}
friedman_to_df <- function(fr, nblocks, nmodels) {
  if (is.null(fr)) return(tibble())
  tibble(
    statistic = unname(fr$statistic),
    df        = unname(fr$parameter),
    p_value   = fr$p.value,
    method    = fr$method,
    n_blocks  = nblocks,
    n_models  = nmodels
  )
}

# --- Overall (pooled) by model from all predictions
cal_overall <- GPPmax_predictions %>%
  group_by(Model) %>%
  summarise(
    n_points           = sum(is.finite(GPPmax_obs) & is.finite(GPPmax_mod)),
    RMSE_cal_overall   = hydroGOF::rmse(GPPmax_obs, GPPmax_mod),
    R2_cal_overall     = adj_r2_from_vectors(GPPmax_obs, GPPmax_mod),
    .groups = "drop"
  )

# AIC/BIC summaries from all fits
aicbic_by_model <- results_all %>%
  group_by(Model) %>%
  summarise(
    n_fits   = n(),
    AIC_mean = mean(AIC, na.rm = TRUE),
    AIC_sum  = sum(AIC,  na.rm = TRUE),
    BIC_mean = mean(BIC, na.rm = TRUE),
    BIC_sum  = sum(BIC,  na.rm = TRUE),
    .groups = "drop"
  )

calibration_summary <- cal_overall %>% left_join(aicbic_by_model, by = "Model")

# --- Blocked comparisons (subject = vegetation group)
perf_subject <- GPPmax_predictions %>%
  group_by(Veg, Model) %>%
  summarise(
    RMSE_subj = hydroGOF::rmse(GPPmax_obs, GPPmax_mod),
    R2_subj   = adj_r2_from_vectors(GPPmax_obs, GPPmax_mod),
    .groups = "drop"
  )

# Friedman + Holm post-hoc (RMSE)
rmse_wide <- perf_subject %>%
  select(Veg, Model, RMSE_subj) %>%
  tidyr::pivot_wider(names_from = Model, values_from = RMSE_subj) %>% na.omit()
friedman_rmse <- NULL; posthoc_rmse <- NULL
if (nrow(rmse_wide) > 0 && (ncol(rmse_wide) - 1) >= 2) {
  rmse_mat <- as.matrix(rmse_wide[, setdiff(names(rmse_wide), "Veg")])
  friedman_rmse <- friedman.test(rmse_mat)
  keep_idx <- perf_subject$Veg %in% rmse_wide$Veg
  posthoc_rmse <- pairwise.wilcox.test(
    x = perf_subject$RMSE_subj[keep_idx],
    g = perf_subject$Model[keep_idx],
    paired = TRUE, p.adjust.method = "holm"
  )
}

# Friedman + Holm post-hoc (R²)
r2_wide <- perf_subject %>%
  select(Veg, Model, R2_subj) %>%
  tidyr::pivot_wider(names_from = Model, values_from = R2_subj) %>% na.omit()
friedman_r2 <- NULL; posthoc_r2 <- NULL
if (nrow(r2_wide) > 0 && (ncol(r2_wide) - 1) >= 2) {
  r2_mat <- as.matrix(r2_wide[, setdiff(names(r2_wide), "Veg")])
  friedman_r2 <- friedman.test(r2_mat)
  keep_idx <- perf_subject$Veg %in% r2_wide$Veg
  posthoc_r2 <- pairwise.wilcox.test(
    x = perf_subject$R2_subj[keep_idx],
    g = perf_subject$Model[keep_idx],
    paired = TRUE, p.adjust.method = "holm"
  )
}

# --- Save all stats to Excel
wb <- createWorkbook()
addWorksheet(wb, "Calibration_overall"); writeData(wb, "Calibration_overall", calibration_summary)
addWorksheet(wb, "Perf_by_Veg");         writeData(wb, "Perf_by_Veg", perf_subject)

if (nrow(rmse_wide) > 0) { addWorksheet(wb, "RMSE_Wide"); writeData(wb, "RMSE_Wide", rmse_wide) }
if (nrow(r2_wide)   > 0) { addWorksheet(wb, "R2_Wide");   writeData(wb, "R2_Wide",   r2_wide)   }

if (!is.null(friedman_rmse)) {
  addWorksheet(wb, "Friedman_RMSE")
  writeData(wb, "Friedman_RMSE", friedman_to_df(friedman_rmse, nrow(rmse_wide), ncol(rmse_wide)-1))
}
if (!is.null(friedman_r2)) {
  addWorksheet(wb, "Friedman_R2")
  writeData(wb, "Friedman_R2", friedman_to_df(friedman_r2, nrow(r2_wide), ncol(r2_wide)-1))
}
if (!is.null(posthoc_rmse)) {
  addWorksheet(wb, "PostHoc_RMSE_Long")
  writeData(wb, "PostHoc_RMSE_Long", matrix_to_long(posthoc_rmse$p.value))
}
if (!is.null(posthoc_r2)) {
  addWorksheet(wb, "PostHoc_R2_Long")
  writeData(wb, "PostHoc_R2_Long", matrix_to_long(posthoc_r2$p.value))
}

saveWorkbook(wb, "GPPmax_Model_Comparison_Stats.xlsx", overwrite = TRUE)
cat("✅ Saved: GPPmax_Model_Comparison_Stats.xlsx (overall metrics + Friedmans + post-hoc)\n")

# =============================
# Observed vs Modeled GPP — Time Series
# Replicate pages (points) + Veg means pages (mean ± SE)
# =============================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

# ---- Define vegetation levels
veg_levels <- c("B","BG","BGS","GS","G")

# ---- Load data with predictions (or use in-memory object)
if (!exists("cbd_with_params")) {
  cbd_with_params <- read.csv("mesocosm_data1_with_GPPmax_predictions.csv",
                              stringsAsFactors = FALSE)
}

# ---- Prepare modeled GPP using the rectangular-hyperbola
# ---- Prepare modeled GPP using the rectangular-hyperbola
# Use replicate-specific K values (provided by user)
k_values <- c(
  "B1" = 407.0161559, "B2" = 396.1684879, "B3" = 744.6314954,
  "BG1" = 1452.782738, "BG2" = 426.3371158, "BG3" = 516.7650617,
  "BGS1" = 342.1892956, "BGS2" = 556.3127428, "BGS3" = 448.493284,
  "G1" = 93.29724454,  "G2" = 621.400516,  "G3" = 959.9439042,
  "GS1" = 381.1992513, "GS2" = 531.9358794, "GS3" = 1059.06864
)

# ---- Compute modeled GPP using replicate-level K values
gpp_df_src <- cbd_with_params %>%
  mutate(
    Plot = as.character(Plot),
    # Match K based on replicate Plot name (e.g., "B1", "BG2", ...)
    K_value = as.numeric(k_values[Plot]),
    PAR  = as.numeric(PAR),
    GPP  = as.numeric(GPP),
    GPPmax_LAI_pred = as.numeric(GPPmax_LAI_pred)
  ) %>%
  mutate(
    # Calculate modeled GPP from GPPmax and PAR using K specific to replicate
    GPP_LAI_modeled = ifelse(
      is.finite(GPPmax_LAI_pred) &
        is.finite(PAR) & is.finite(K_value),
      (GPPmax_LAI_pred * PAR) / (K_value + PAR),
      NA_real_
    )
  )


# ---- Replicate-level long table (Observed / Modeled)
ts_long_repl <- gpp_df_src %>%
  select(Date, Plot, GPP, GPP_LAI_modeled) %>%
  pivot_longer(cols = c(GPP, GPP_LAI_modeled),
               names_to = "Type", values_to = "Value") %>%
  mutate(
    Type = recode(Type, GPP = "Observed", GPP_LAI_modeled = "Modeled"),
    Type = factor(Type, levels = c("Observed", "Modeled"))
  )

if (any(is.na(gpp_df_src$K_value))) {
  missing_plots <- unique(gpp_df_src$Plot[is.na(gpp_df_src$K_value)])
  warning("⚠️ Missing K values for plots: ", paste(missing_plots, collapse = ", "))
} else {
  cat("✅ All replicate K values matched successfully.\n")
}

# ---- Vegetation-mean ± SE across replicates
se_fun <- function(x) {
  n <- sum(is.finite(x))
  if (n <= 1) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / sqrt(n)
}

veg_agg <- gpp_df_src %>%
  group_by(Plot, Date) %>%
  summarise(
    mean_Obs = mean(GPP, na.rm = TRUE),
    se_Obs   = se_fun(GPP),
    mean_Mod = mean(GPP_LAI_modeled, na.rm = TRUE),
    se_Mod   = se_fun(GPP_LAI_modeled),
    .groups  = "drop"
  ) %>%
  pivot_longer(cols = c(mean_Obs, se_Obs, mean_Mod, se_Mod),
               names_to = c(".value","Type"),
               names_pattern = "(mean|se)_(Obs|Mod)") %>%
  mutate(
    Type  = factor(recode(Type, Obs = "Observed", Mod = "Modeled"),
                   levels = c("Observed", "Modeled"))
  )

# ---- Plot helpers
plot_ts_replicates <- function(df, title_prefix) {
  ggplot(df, aes(x = as.POSIXct(Date), y = Value, color = Type)) +
    geom_point(size = 1.5, alpha = 0.9, na.rm = TRUE) +
    facet_wrap(~ Plot, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
    labs(
      title = sprintf("%s — Replicates (Observed vs Modeled)", title_prefix),
      x = "Date",
      y = expression("GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
      color = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "top")
}

plot_ts_means <- function(df, title_prefix) {
  ggplot(df, aes(x = as.POSIXct(Date), y = mean, color = Type)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0, alpha = 0.7, na.rm = TRUE) +
    geom_point(size = 1.8, alpha = 0.9, na.rm = TRUE) +
    facet_wrap(~ Plot, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
    labs(
      title = sprintf("%s — Vegetation Means ± SE", title_prefix),
      x = "Date",
      y = expression("GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
      color = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "top")
}

# ---- Create the PDF (2 pages)
pdf("Obs_vs_Modeled_GPP_TimeSeries_GroupLevel.pdf", width = 14, height = 9)
print(plot_ts_replicates(ts_long_repl, "Observed vs Modeled GPP"))
print(plot_ts_means(veg_agg, "Observed vs Modeled GPP"))
dev.off()
cat("✅ Saved: Obs_vs_Modeled_GPP_TimeSeries_GroupLevel.pdf (replicates & veg means)\n")

# =============================
#   GPPmax time-series: Observed vs Modeled (POINTS ONLY)
# =============================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

gppmax_ts <- GPPmax_predictions %>%
  mutate(
    Date   = as.POSIXct(Date, format = "%m/%d/%Y %H:%M"),
    Veg    = factor(Veg, levels = group_levels_order, ordered = TRUE),
    Model  = factor(Model, levels = c("Model_LAI", "Model_VI")),
    Label  = ifelse(Model == "Model_LAI", "LAI", "VI")
  ) %>%
  filter(Label == "LAI" | (USE_VI & Label == "VI"))

# replicate-level long table (Observed / Modeled)
gppmax_repl_long <- gppmax_ts %>%
  select(Date, Veg, Model, Label, GPPmax_obs, GPPmax_mod) %>%
  pivot_longer(cols = c(GPPmax_obs, GPPmax_mod),
               names_to = "Type", values_to = "GPPmax") %>%
  mutate(
    Type = recode(Type, GPPmax_obs = "Observed", GPPmax_mod = "Modeled"),
    Type = factor(Type, levels = c("Observed", "Modeled"))
  )

# vegetation mean per timestamp (average across replicates if multiple entries exist)
gppmax_veg_long <- gppmax_ts %>%
  group_by(Veg, Date, Label) %>%
  summarise(
    Obs = mean(GPPmax_obs, na.rm = TRUE),
    Mod = mean(GPPmax_mod, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Obs, Mod), names_to = "Type", values_to = "GPPmax") %>%
  mutate(
    Type = recode(Type, Obs = "Observed", Mod = "Modeled"),
    Type = factor(Type, levels = c("Observed", "Modeled"))
  )

# ---- Plot helpers — POINTS ONLY ----
plot_gppmax_repl_pts <- function(df_ts, model_label, title_prefix) {
  ggplot(dplyr::filter(df_ts, Label == model_label),
         aes(x = Date, y = GPPmax, color = Type)) +
    geom_point(na.rm = TRUE, size = 1.6, alpha = 0.9) +
    facet_wrap(~ Veg, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(Observed = "black", Modeled = "red")) +
    labs(
      title = sprintf("%s — %s (GPPmax; points per vegetation group)", title_prefix, model_label),
      x = "Date",
      y = expression("GPPmax ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
      color = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "top")
}

plot_gppmax_veg_pts <- function(df_ts, model_label, title_prefix) {
  ggplot(dplyr::filter(df_ts, Label == model_label),
         aes(x = Date, y = GPPmax, color = Type)) +
    geom_point(na.rm = TRUE, size = 1.8, alpha = 0.9) +
    facet_wrap(~ Veg, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(Observed = "black", Modeled = "red")) +
    labs(
      title = sprintf("%s — %s (GPPmax; vegetation means, points)", title_prefix, model_label),
      x = "Date",
      y = expression("GPPmax ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
      color = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "top")
}

# ---- Build the PDF(s)
if (USE_VI) {
  pdf("GPPmax_Obs_vs_Mod_TimeSeries_POINTS_VI_LAI.pdf", width = 14, height = 9)
  print(plot_gppmax_repl_pts(gppmax_repl_long, "LAI", "Observed vs Modeled GPPmax"))  # Page 1
  print(plot_gppmax_repl_pts(gppmax_repl_long, "VI",  "Observed vs Modeled GPPmax"))  # Page 2
  print(plot_gppmax_veg_pts(gppmax_veg_long,  "LAI", "Observed vs Modeled GPPmax"))   # Page 3
  print(plot_gppmax_veg_pts(gppmax_veg_long,  "VI",  "Observed vs Modeled GPPmax"))   # Page 4
  dev.off()
  message("✅ Saved: GPPmax_Obs_vs_Mod_TimeSeries_POINTS_VI_LAI.pdf")
} else {
  pdf("GPPmax_Obs_vs_Mod_TimeSeries_POINTS_LAI_only.pdf", width = 14, height = 9)
  print(plot_gppmax_repl_pts(gppmax_repl_long, "LAI", "Observed vs Modeled GPPmax"))  # Page 1
  print(plot_gppmax_veg_pts(gppmax_veg_long,  "LAI", "Observed vs Modeled GPPmax"))   # Page 2
  dev.off()
  message("✅ Saved: GPPmax_Obs_vs_Mod_TimeSeries_POINTS_LAI_only.pdf")
}

