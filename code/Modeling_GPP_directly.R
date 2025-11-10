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
library(broom)
library(gridExtra)
library(ggforce)
library(patchwork)
library(stringr)

# =============================
#        Load and Prepare Data
# =============================
cbd <- read_excel("C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by 2nd method/by treatment/other model/mesocosm_data_GPP.xlsx")
cbd <- dplyr::rename(cbd, Date = TIMESTAMP.Reco)
cbd$Date <- as.POSIXct(cbd$Date, format = "%m/%d/%Y")

# Split calibration and evaluation
GPP_cal  <- subset(cbd, ID_camp %in% c("1","3","4","6","8","10","11","13","15","16","18","19","21","22","24","25","27","28","30"))
GPP_eval <- subset(cbd, ID_camp %in% c("2","5","7","9","12","14","17","20","23","26"))

  # Kelvin

setwd("C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by 2nd method/by treatment/other model/for git")

Tmax = 40
Tmin = 0
# =============================
#   CONFIG
# =============================
veg_types <- c("B","BG","BGS","GS","G")
num_iter  <- 1

results_all      <- data.frame()  # calibration fit stats + params
fitted_models    <- list()        # nlsLM objects
params_list      <- list()        # coefficients by model
GPP_predictions <- data.frame()  # calibration predictions (obs vs mod)

# =============================
#   Manual prediction helper (LAI ONLY; NO clamping)
# =============================

predict_manual <- function(fit, newdata) {
  cf   <- coef(fit)
  a    <- cf["a"]
  b    <- cf["b"]
  Topt <- cf["Topt"]
  k    <- cf["k"]   # <--- uppercase, same as model
  
  Ta   <- newdata$Ta
  LAI  <- newdata$LAI
  PAR  <- newdata$PAR
  
  num  <- (Ta - Tmin) * (Ta - Tmax)
  den  <- num - (Ta - Topt)^2
  temp_resp <- ifelse(is.finite(den) & den != 0, num / den, 0)
  
  base_term <- (a * LAI) + b
  
  (base_term * temp_resp * PAR) / (k + PAR)
}



needed_vars_for_model <- function() c("GPP","Ta","LAI", "PAR")

fit_models_for_veg <- function(veg_name, GPP_cal, num_iter = 1) {
  veg_data <- subset(GPP_cal, Plot == veg_name)  # all reps of this plot
  
  for (i in 1:num_iter) {
    tryCatch({
      need_vars <- needed_vars_for_model()
      
      rep_data_f <- veg_data %>%
        mutate(across(all_of(need_vars), as.numeric)) %>%
        filter(if_all(all_of(need_vars), ~ !is.na(.)))
      
      if (nrow(rep_data_f) == 0) stop("No valid rows for fitting")
      
      # ===============================
      # model (same as yours, just tidy)
      # ===============================
      formula_model <- GPP ~ (
        ((a * LAI) + b) *
          (((Ta - Tmin) * (Ta - Tmax)) /
             (((Ta - Tmin) * (Ta - Tmax)) - (Ta - Topt)^2)) *
          PAR
      ) / (k + PAR)
      
      # ------------------------------------------------
      # IMPORTANT: give starts for *all* free parameters
      # ------------------------------------------------
      # tweak these if convergence is slow
      start_list <- list(
        a    = 0.5,
        b    = 0.1,   # °C, upper bound of response
        Topt = 20,    # °C, optimum
        k    = 100    # µmol m-2 s-1, half-sat for PAR part
      )
      
      # you can keep them wide open:
      lower_bounds <- c(
        a    = -Inf,
        b    = -Inf, # must be > Tmin later, but nlsLM can handle it
        Topt = 0,
        k    = 10
      )
      upper_bounds <- c(
        a    =  Inf,
        b    =  Inf,
        Topt = 40,
        k    = 2000
      )
      
      fit <- nlsLM(
        formula_model,
        data   = rep_data_f,
        start  = start_list,
        lower  = lower_bounds,
        upper  = upper_bounds,
        control = nls.lm.control(maxiter = 1000, ftol = 1e-10)
      )
      
      # ---------------------------------
      # stats
      # ---------------------------------
      preds    <- predict(fit)
      adj_r2   <- summary(lm(rep_data_f$GPP ~ preds))$adj.r.squared
      rmse_val <- rmse(rep_data_f$GPP, preds)
      aic_val  <- AIC(fit)
      bic_val  <- BIC(fit)
      
      model_key <- paste0(veg_name, "_iter", i, "_Model_LAI")
      fitted_models[[model_key]] <<- fit
      params_list[[model_key]]  <<- coef(fit)
      
      sum_fit <- summary(fit)
      coefs   <- coef(sum_fit)
      params  <- as.numeric(coefs[, "Estimate"])
      errors  <- as.numeric(coefs[, "Std. Error"])
      t_stat  <- params / errors
      df_fit  <- df.residual(fit)
      p_vals  <- 2 * pt(-abs(t_stat), df = df_fit)
      
      row <- list(
        SubsetType = veg_name,
        Iteration  = i,
        Model      = "Model_LAI",
        AIC        = aic_val,
        BIC        = bic_val,
        RMSE       = rmse_val,
        R2         = adj_r2,
        Veg        = veg_name
      )
      
      for (j in seq_along(params)) {
        pname <- rownames(coefs)[j]
        row[[pname]]                   <- params[j]
        row[[paste0("SD_", pname)]]    <- errors[j]
        row[[paste0("pval_", pname)]]  <- p_vals[j]
      }
      
      results_all <<- bind_rows(results_all, as.data.frame(row))
      
      GPP_predictions <<- bind_rows(
        GPP_predictions,
        data.frame(
          Veg       = veg_name,
          Rep       = NA_character_,
          Model     = "Model_LAI",
          Date      = rep_data_f$Date,
          GPP_obs   = rep_data_f$GPP,
          GPP_mod   = preds,
          Residual  = rep_data_f$GPP - preds
        )
      )
      
    }, error = function(e) {
      cat("❌ Error in", veg_name, "Model_LAI :", e$message, "\n")
    })
  }
}

# =============================
#   Run fitting (LAI only, by Plot)
# =============================
for (veg in veg_types) fit_models_for_veg(veg, GPP_cal, num_iter)

# =============================
#   SAVE results_all EXPANDED TO REPLICATES (B1..G3)
# =============================
R_gas <- 8.314   # J mol^-1 K^-1
TQ    <- 288.15  # K, for Q10

results_all_repl <- results_all %>%
  # duplicate each row for Rep = 1,2,3 and stamp SubsetType as B1..G3
  tidyr::crossing(Rep = c("1","2","3")) %>%
  dplyr::mutate(
    SubsetType = paste0(Veg, Rep)
  ) %>%
  # keep column order sensible: SubsetType then everything else
  dplyr::relocate(SubsetType, .before = dplyr::everything())

openxlsx::write.xlsx(results_all_repl, "GPP_Model_Results_All.xlsx", rowNames = FALSE)
cat("✅ GPP_Model_Results_All.xlsx saved with parameters duplicated for B1..G3, BG1..BG3, ..., G1..G3.\n")


# =============================
#   Calibration summary table (with overall p-value)
# =============================
cal_results <- GPP_predictions %>%
  group_by(Veg, Model) %>%
  summarise(
    {
      keep <- is.finite(GPP_obs) & is.finite(GPP_mod)
      if (sum(keep) < 3) {
        tibble::tibble(
          RMSE_cal = NA_real_,
          R2_cal   = NA_real_,
          P_cal    = NA_real_
        )
      } else {
        y_obs <- GPP_obs[keep]; y_mod <- GPP_mod[keep]
        lin   <- lm(y_obs ~ y_mod)
        tibble::tibble(
          RMSE_cal = hydroGOF::rmse(y_obs, y_mod),
          R2_cal   = summary(lin)$adj.r.squared,
          P_cal    = broom::glance(lin)$p.value   # overall model p-value (F-test)
        )
      }
    },
    .groups = "drop"
  ) %>%
  mutate(Rep = NA_character_) %>%
  relocate(Rep, .after = Veg)

openxlsx::write.xlsx(cal_results, "GPP_Model_Calibration_Metrics.xlsx", rowNames = FALSE)

# =============================
#   Evaluation (LAI only, by Plot) — with overall p-value
# =============================
eval_results <- data.frame()

for (model_key in names(fitted_models)) {
  fit <- fitted_models[[model_key]]
  mi <- strcapture("^([A-Z]+)_iter[0-9]+_(.+)$", model_key,
                   proto = list(Veg = character(), Model = character()))
  veg_i <- mi$Veg
  
  need <- needed_vars_for_model()
  eval_ready <- subset(GPP_eval, Plot == veg_i) %>%
    mutate(across(all_of(need), as.numeric)) %>%
    filter(if_all(all_of(need), ~ is.finite(.)))
  if (nrow(eval_ready) == 0) { cat("\u26A0\ufe0f No eval data for", model_key, "\n"); next }
  
  preds_eval <- predict_manual(fit, eval_ready)
  keep <- is.finite(eval_ready$GPP) & is.finite(preds_eval)
  if (sum(keep) < 3) { cat("\u26A0\ufe0f Insufficient valid points for", model_key, "\n"); next }
  
  y_obs <- eval_ready$GPP[keep]; y_mod <- preds_eval[keep]
  lin   <- lm(y_obs ~ y_mod)
  
  eval_results <- rbind(
    eval_results,
    data.frame(
      Veg        = veg_i,
      Rep        = NA_character_,
      Model      = "Model_LAI",
      ModelKey   = model_key,
      R2_eval    = summary(lin)$adj.r.squared,
      RMSE_eval  = hydroGOF::rmse(y_obs, y_mod),
      P_eval     = broom::glance(lin)$p.value   # overall model p-value (F-test)
    )
  )
}

openxlsx::write.xlsx(eval_results, "GPP_Model_Evaluation_Metrics.xlsx", rowNames = FALSE)

# =============================
#   Join calibration & evaluation, rank (LAI only)
# =============================
cal_tidy <- results_all %>%
  mutate(Rep = NA_character_) %>%
  select(Veg, Rep, Model, RMSE_cal = RMSE, R2_cal = R2, AIC, BIC)

perf_all <- eval_results %>% left_join(cal_tidy, by = c("Veg","Rep","Model"))

perf_ranked <- perf_all %>%
  group_by(Veg) %>%
  arrange(RMSE_eval, desc(R2_eval), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

openxlsx::write.xlsx(perf_all,    "GPP_Model_Performance_Cal_vs_Eval_LAI.xlsx", rowNames = FALSE)
openxlsx::write.xlsx(perf_ranked, "GPP_Model_Performance_Ranked_LAI.xlsx",      rowNames = FALSE)

# =============================
#   Build Evaluation predictions for saving
# =============================
GPP_predictions$Set <- "Calibration"

eval_preds_df <- data.frame()
for (model_key in names(fitted_models)) {
  fit <- fitted_models[[model_key]]
  mi <- strcapture("^([A-Z]+)_iter[0-9]+_(.+)$",
                   model_key,
                   proto = list(Veg = character(), Model = character()))
  veg_i <- mi$Veg
  
  need <- needed_vars_for_model()
  subset_eval <- GPP_eval %>%
    dplyr::filter(Plot == veg_i) %>%
    dplyr::mutate(across(all_of(need), as.numeric)) %>%
    dplyr::filter(if_all(all_of(need), ~ !is.na(.)))
  if (nrow(subset_eval) == 0) next
  
  preds_eval <- predict_manual(fit, subset_eval)
  
  eval_preds_df <- dplyr::bind_rows(
    eval_preds_df,
    data.frame(
      Veg = veg_i, Rep = NA_character_, Model = "Model_LAI",
      Date = subset_eval$Date,
      GPP_obs = subset_eval$GPP,
      GPP_mod = preds_eval,
      Residual = subset_eval$GPP - preds_eval,
      Set = "Evaluation",
      stringsAsFactors = FALSE
    )
  )
}

# =============================
#   UPDATED: Save Calibration + Evaluation + Combined
# =============================
make_sheet <- function(df) {
  if (nrow(df) == 0) return(df)
  df %>%
    mutate(
      Plot_ID = Veg,                     # now Plot-level
      Date = format(as.POSIXct(Date), "%m/%d/%Y %H:%M")
    ) %>%
    select(Veg, Rep, Plot_ID, Model, Set, Date, GPP_obs, GPP_mod, Residual)
}
cal_to_save  <- make_sheet(GPP_predictions)
eval_to_save <- make_sheet(eval_preds_df)
comb_to_save <- make_sheet(bind_rows(GPP_predictions, eval_preds_df))
openxlsx::write.xlsx(
  x = list(Calibration = cal_to_save, Evaluation = eval_to_save, Combined = comb_to_save),
  file = "GPP_Observed_vs_Modeled_All_LAI.xlsx",
  overwrite = TRUE
)

# =============================
#   Prepare Combined Dataset for scatter panels (LAI only)
# =============================
combined_preds <- bind_rows(GPP_predictions, eval_preds_df)
combined_preds$Plot_ID <- factor(
  combined_preds$Veg,
  levels = c("B","BG","BGS","GS","G")
)

safe_stats <- function(df) {
  ok <- is.finite(df$GPP_obs) & is.finite(df$GPP_mod)
  df <- df[ok, , drop = TRUE]
  if (nrow(df) < 2 || stats::sd(df$GPP_obs) == 0) {
    return(tibble::tibble(R2 = NA_real_, Pval = NA_real_))
  }
  m <- lm(GPP_mod ~ GPP_obs, data = df)
  tibble::tibble(
    R2   = summary(m)$adj.r.squared,
    Pval = broom::glance(m)$p.value
  )
}

stats_df <- combined_preds %>%
  dplyr::group_split(Plot_ID, Set, Model) %>%
  lapply(function(df) {
    out <- safe_stats(df)
    cbind(
      data.frame(
        Plot_ID = unique(df$Plot_ID),
        Set     = unique(df$Set),
        Model   = unique(df$Model)
      ),
      out
    )
  }) %>% dplyr::bind_rows()
openxlsx::write.xlsx(stats_df, "GPP_Panel_Stats_R2_Pval_LAI.xlsx", rowNames = FALSE)

stats_df <- combined_preds %>%
  group_by(Plot_ID, Set, Model) %>%
  summarise(
    R2 = summary(lm(GPP_mod ~ GPP_obs))$adj.r.squared,
    Pval = glance(lm(GPP_mod ~ GPP_obs))$p.value,
    .groups = "drop"
  )
combined_preds <- left_join(combined_preds, stats_df, by = c("Plot_ID","Set","Model"))

plot_scatter_by_set <- function(set_type) {
  df <- combined_preds %>% filter(Set == set_type, Model == "Model_LAI")
  annotation_df <- df %>%
    group_by(Plot_ID) %>% slice(1) %>%
    mutate(label = paste0("R² = ", round(R2, 2), "\nP = ", signif(Pval, 2)))
  ggplot(df, aes(x = GPP_obs, y = GPP_mod)) +
    geom_point(alpha = 0.6, color = "black", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, color = "gray40", linetype = "dashed") +
    geom_text(data = annotation_df,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, size = 4.5, hjust = 1.1, vjust = 1.5) +
    facet_wrap(~ Plot_ID, ncol = 3, scales = "free") +
    labs(title = paste(set_type, "- Model_LAI - Observed vs Modeled GPP (by Plot)"),
         x = expression(Observed~GPP~(mu*mol~m^{-2}~s^{-1})),
         y = expression(Modeled~GPP~(mu*mol~m^{-2}~s^{-1}))) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.text = element_text(size = 9),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank())
}
# Save scatter panels
fig_cal  <- plot_scatter_by_set("Calibration")
fig_eval <- plot_scatter_by_set("Evaluation")
ggsave("GPP_Calibration_Model_LAI_Scatterplots_byReplicate.png", fig_cal,  width = 10, height = 7, dpi = 300)
ggsave("GPP_Evaluation_Model_LAI_Scatterplots_byReplicate.png",  fig_eval, width = 10, height = 7, dpi = 300)

# =============================
#   Time Series Plots (calibration set only) — LAI only
# =============================
data <- GPP_predictions
data$Date <- as.POSIXct(data$Date)
data$Plot_ID <- factor(data$Veg, levels = c("B","BG","BGS","GS","G"))

data_long <- data %>%
  pivot_longer(cols = c(GPP_obs, GPP_mod), names_to = "Type", values_to = "GPP") %>%
  filter(!is.na(GPP)) %>%
  mutate(Type = recode(Type, GPP_obs = "Observed", GPP_mod = "Modeled"))

n_col <- 3; n_row <- 2
plots_per_page <- n_col * n_row
total_panels <- length(unique(data_long$Plot_ID))
n_pages <- ceiling(total_panels / plots_per_page)

for (i in 1:n_pages) {
  p <- ggplot(data_long, aes(x = Date, y = GPP, color = Type)) +
    geom_point(size = 1.5) +
    ggforce::facet_grid_paginate(~ Plot_ID, ncol = n_col, nrow = n_row, page = i, scales = "free_y") +
    labs(title = paste("Observed vs Modeled (Model_LAI) - Page", i),
         x = "Date", y = expression(~(mu*mol~m^{-2}~s^{-1}))) +
    scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(size = 9),
          legend.position = "top")
  print(p)
}
dev.off()
cat("✅ Paginated time series plot saved (Model_LAI).\n")

plot_all_replicates_by_model <- function() {
  df <- data_long
  file_name <- "_Observed_vs_Modeled_AllReplicates_Model_LAI.png"
  png(file_name, width = 12, height = 8, units = "in", res = 300)
  print(
    ggplot(df, aes(x = Date, y = GPP, color = Type)) +
      geom_point(size = 1.4) +
      facet_wrap(~ Plot_ID, ncol = 3, scales = "free_y") +
      labs(title = "Observed vs Modeled – All Plots (Model_LAI)",
           x = "Date", y = expression(~(mu*mol~m^{-2}~s^{-1}))) +
      scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
      theme_bw(base_size = 14) +
      theme(strip.text = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 90, size = 9),
            axis.text.y = element_text(size = 9),
            legend.position = "top")
  )
  dev.off(); cat("✅ Saved:", file_name, "\n")
}
plot_all_replicates_by_model()

# =============================
#   SIGNIFICANCE TESTS & OVERALL METRICS (LAI only)
# =============================
adj_r2_from_vectors <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (!any(keep)) return(NA_real_)
  summary(lm(y[keep] ~ yhat[keep]))$adj.r.squared
}
matrix_to_long <- function(mat) {
  if (is.null(mat)) return(tibble(Model1=character(), Model2=character(), p_adj=numeric()))
  as_tibble(mat, rownames = "Model1") %>%
    pivot_longer(-Model1, names_to = "Model2", values_to = "p_adj") %>%
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

cal_overall <- GPP_predictions %>%
  group_by(Model) %>%
  summarise(
    n_points         = sum(is.finite(GPP_obs) & is.finite(GPP_mod)),
    RMSE_cal_overall = hydroGOF::rmse(GPP_obs, GPP_mod),
    R2_cal_overall   = adj_r2_from_vectors(GPP_obs, GPP_mod),
    .groups = "drop"
  )

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

eval_summary <- eval_results %>%
  group_by(Model) %>%
  summarise(
    n_reps           = n(),
    RMSE_eval_mean   = mean(RMSE_eval, na.rm = TRUE),
    RMSE_eval_median = median(RMSE_eval, na.rm = TRUE),
    R2_eval_mean     = mean(R2_eval, na.rm = TRUE),
    R2_eval_median   = median(R2_eval, na.rm = TRUE),
    .groups = "drop"
  )

wb <- createWorkbook()
addWorksheet(wb, "Calibration_overall"); writeData(wb, "Calibration_overall", calibration_summary)
addWorksheet(wb, "Evaluation_overall");  writeData(wb, "Evaluation_overall",  eval_summary)
saveWorkbook(wb, "Model_Comparison_Stats_LAI_only.xlsx", overwrite = TRUE)
cat("✅ Saved: Model_Comparison_Stats_LAI_only.xlsx (overall metrics)\n")

# ===================================================================
#  PARAMETER PLOTS — LAI only (Ea_kJmol with SE & stars)
# ===================================================================
R_gas <- 8.314  # J mol^-1 K^-1
params <- readxl::read_excel("GPP_Model_Results_All.xlsx")

params <- params %>%
  mutate(
    SubsetType = Veg                    # ensure SubsetType is the Plot label
  )

subset_levels <- c("B","BG","BGS","GS","G")

params_tidy <- params %>%
  transmute(
    SubsetType, Model,
    a, b, Topt, k,
    SD_a = dplyr::coalesce(SD_a, NA_real_),
    SD_b = dplyr::coalesce(SD_b, NA_real_),
    SD_Topt = dplyr::coalesce(SD_Topt, NA_real_),
    SD_k = dplyr::coalesce(SD_k, NA_real_),
    pval_a = dplyr::coalesce(pval_a, NA_real_),
    pval_b = dplyr::coalesce(pval_b, NA_real_),
    pval_Topt = dplyr::coalesce(pval_Topt, NA_real_),
    pval_k = dplyr::coalesce(pval_k, NA_real_)
  ) %>%
  pivot_longer(cols = c(a, b, Topt, k), names_to = "Parameter", values_to = "Estimate") %>%
  mutate(
    SE = dplyr::case_when(
      Parameter == "a"        ~ SD_a,
      Parameter == "b"        ~ SD_b,
      Parameter == "Topt"        ~ SD_Topt,
      Parameter == "k" ~ SD_k
    ),
    Pval = dplyr::case_when(
      Parameter == "a"        ~ pval_a,
      Parameter == "b"        ~ pval_b,
      Parameter == "Topt"        ~ pval_Topt,
      Parameter == "k" ~ pval_k
    ),
    Signif     = ifelse(Pval < 0.05, "*", ""),
    SubsetType = factor(SubsetType, levels = subset_levels),
    Model      = factor("Model_LAI", levels = "Model_LAI")
  )

plot_param_hist <- function(param_name, ylab_default = "Estimate ± SE") {
  df <- dplyr::filter(params_tidy, Parameter == param_name)
  if (nrow(df) == 0) { warning(paste("No data found for", param_name)); return(NULL) }
  ylab <- if (param_name == "Ea_kJmol") "Estimate ± SE (kJ mol^-1)" else ylab_default
  ggplot(df, aes(x = SubsetType, y = Estimate)) +
    geom_col(color = "black", fill = NA, width = 0.7) +
    geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2) +
    geom_text(aes(label = Signif, y = Estimate + SE + 0.1), vjust = 0, size = 5) +
    facet_grid(cols = vars(Model), scales = "free_y") +
    labs(title = paste("Parameter", param_name, "by Plot (Model_LAI)"),
         y = ylab, x = "Plot") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 16),
          strip.text = element_text(size = 12, face = "bold"),
          legend.position = "none")
}

pdf("Parameter_a_b_c_Ea_kJmol_Histograms_LAI.pdf", width = 14, height = 8)
for (p in c("a","b","Topt","k")) {
  p_obj <- plot_param_hist(p)
  if (!is.null(p_obj)) print(p_obj)
}
dev.off()
cat("✅ Per-plot parameter histograms saved (Model_LAI).\n")

# =============================
#  Two-page PDFs (Cal & Eval) — LAI only
#  Pages:
#   1) Plot-level (B,BG,BGS,GS,G)
#   2) Same as 1 (kept structure)
# =============================

build_set_df <- function(df) {
  req <- c("Veg","Model","Date","GPP_obs","GPP_mod")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing columns in input df: ", paste(miss, collapse = ", "))
  df %>%
    mutate(
      Date   = as.POSIXct(Date),
      Veg    = as.character(Veg),
      Rep    = NA_character_,
      Model  = "Model_LAI",
      Plot_ID = Veg
    )
}

cal_df  <- build_set_df(GPP_predictions)
eval_df <- build_set_df(eval_preds_df)

rep_levels <- c("B","BG","BGS","GS","G")

plot_replicates_3x5 <- function(df, title_prefix) {
  df_m <- df %>%
    filter(Model == "Model_LAI") %>%
    mutate(Plot_ID = factor(Plot_ID, levels = rep_levels)) %>%
    filter(Plot_ID %in% rep_levels)
  if (nrow(df_m) == 0) return(ggplot() + ggtitle(paste(title_prefix, "(no data)")) + theme_void())
  df_long <- df_m %>%
    pivot_longer(c(GPP_obs, GPP_mod), names_to = "Type", values_to = "GPP") %>%
    mutate(Type = recode(Type, GPP_obs = "Observed", GPP_mod = "Modeled"))
  ggplot(df_long, aes(x = Date, y = GPP, color = Type)) +
    geom_point(size = 1.2, na.rm = TRUE) +
    facet_wrap(~ Plot_ID, ncol = 3, nrow = 2, scales = "free_y") +
    labs(
      title = paste(title_prefix, "— Model_LAI (Plots)"),
      x = "Date", y = expression(~(mu*mol~m^{-2}~s^{-1}))
    ) +
    scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
      legend.position = "top"
    )
}

plot_group_avg_3x2 <- function(df, title_prefix) {
  # Already at plot-level; reuse same layout
  plot_replicates_3x5(df, paste0(title_prefix, " (avg layout)"))
}

make_two_page_pdf <- function(df, file_out, title_prefix) {
  p1 <- plot_replicates_3x5(df, title_prefix)  # Page 1
  p2 <- plot_group_avg_3x2(df, title_prefix)   # Page 2
  pdf(file_out, width = 14, height = 10); print(p1); print(p2); dev.off()
  message("✅ Saved: ", file_out)
}

make_two_page_pdf(cal_df,  "_Observed_vs_Modeled_CALIBRATION_Model_LAI_2pages.pdf", "Calibration")
make_two_page_pdf(eval_df, "_Observed_vs_Modeled_EVALUATION_Model_LAI_2pages.pdf",  "Evaluation")

# ============================================================
# EXTRA: Replicate-view outputs using plot-level fitted models
# Files created (replicate aspect):
#   _Observed_vs_Modeled_AllReplicates_Model_LAI_BY_REP.png
#   _Observed_vs_Modeled_CALIBRATION_Model_LAI_2pages_BY_REP.pdf
#   _Observed_vs_Modeled_EVALUATION_Model_LAI_2pages_BY_REP.pdf
#   GPP_Calibration_Model_LAI_Scatterplots_byReplicate_BY_REP.png
#   GPP_Evaluation_Model_LAI_Scatterplots_byReplicate_BY_REP.png
#   RMSE_vs_R2_Calibration_and_Evaluation_Model_LAI_BY_REP.png
# ============================================================

# Helper to extract Veg from model_key built earlier (B, BG, BGS, GS, G)
.parse_key_to_veg <- function(model_key) {
  mi <- strcapture("^([A-Z]+)_iter[0-9]+_(.+)$", model_key,
                   proto = list(Veg = character(), Model = character()))
  mi$Veg
}

rep_levels_15 <- c("B1","B2","B3",
                   "BG1","BG2","BG3",
                   "BGS1","BGS2","BGS3",
                   "GS1","GS2","GS3",
                   "G1","G2","G3")

# ---------- Build per-replicate Calibration & Evaluation predictions ----------
build_replicate_preds <- function(base_df, set_name) {
  out <- list()
  for (model_key in names(fitted_models)) {
    fit  <- fitted_models[[model_key]]
    veg  <- .parse_key_to_veg(model_key)
    need <- needed_vars_for_model()
    for (rep_i in c("1","2","3")) {
      df <- base_df %>%
        dplyr::filter(Plot == veg, Rep == rep_i) %>%
        dplyr::mutate(across(all_of(need), as.numeric)) %>%
        dplyr::filter(if_all(all_of(need), ~ is.finite(.)))
      if (!nrow(df)) next
      preds <- predict_manual(fit, df)
      out[[length(out) + 1]] <- data.frame(
        Veg      = veg,
        Rep      = rep_i,
        Plot_ID  = paste0(veg, rep_i),
        Model    = "Model_LAI",
        Date     = df$Date,
        GPP_obs = df$GPP,
        GPP_mod = preds,
        Residual = df$GPP - preds,
        Set      = set_name,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out)) dplyr::bind_rows(out) else
    data.frame(Veg=character(), Rep=character(), Plot_ID=character(),
               Model=character(), Date=as.POSIXct(character()),
               GPP_obs=numeric(), GPP_mod=numeric(), Residual=numeric(),
               Set=character(), stringsAsFactors=FALSE)
}

cal_byrep  <- build_replicate_preds(GPP_cal,  "Calibration")
eval_byrep <- build_replicate_preds(GPP_eval, "Evaluation")

# ---------- SCATTERPLOTS (BY REPLICATE) ----------
scatter_byrep <- function(df, title_txt) {
  # stats per panel
  stats_df <- df %>%
    group_by(Plot_ID) %>%
    summarise(
      R2   = { m <- try(lm(GPP_mod ~ GPP_obs, data = .), silent = TRUE);
      if (inherits(m, "try-error")) NA_real_ else summary(m)$adj.r.squared },
      Pval = { m <- try(lm(GPP_mod ~ GPP_obs, data = .), silent = TRUE);
      if (inherits(m, "try-error")) NA_real_ else broom::glance(m)$p.value },
      .groups = "drop"
    )
  df2 <- df %>% mutate(Plot_ID = factor(Plot_ID, levels = rep_levels_15))
  ann <- stats_df %>%
    mutate(label = paste0("R² = ", round(R2, 2), "\nP = ", signif(Pval, 2)))
  ggplot(df2, aes(x = GPP_obs, y = GPP_mod)) +
    geom_point(alpha = 0.6, size = 1.2, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 4.2) +
    facet_wrap(~ Plot_ID, ncol = 3, scales = "free") +
    labs(title = title_txt,
         x = expression(Observed~GPP~(mu*mol~m^{-2}~s^{-1})),
         y = expression(Modeled~GPP~(mu*mol~m^{-2}~s^{-1}))) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank())
}

if (nrow(cal_byrep)) {
  p_cal_rep <- scatter_byrep(cal_byrep,  "Calibration — Model_LAI — Observed vs Modeled (by Replicate)")
  ggsave("GPP_Calibration_Model_LAI_Scatterplots_byReplicate_BY_REP.png",
         p_cal_rep, width = 10, height = 15, dpi = 300)
}
if (nrow(eval_byrep)) {
  p_eval_rep <- scatter_byrep(eval_byrep, "Evaluation — Model_LAI — Observed vs Modeled (by Replicate)")
  ggsave("GPP_Evaluation_Model_LAI_Scatterplots_byReplicate_BY_REP.png",
         p_eval_rep, width = 10, height = 15, dpi = 300)
}

# ---------- TIME SERIES PANELS (BY REPLICATE) ----------
timeseries_byrep <- function(df, title_prefix) {
  df2 <- df %>%
    mutate(Plot_ID = factor(Plot_ID, levels = rep_levels_15)) %>%
    tidyr::pivot_longer(c(GPP_obs, GPP_mod), names_to = "Type", values_to = "GPP") %>%
    mutate(Type = recode(Type, GPP_obs = "Observed", GPP_mod = "Modeled"))
  ggplot(df2, aes(x = Date, y = GPP, color = Type)) +
    geom_point(size = 1.2, na.rm = TRUE) +
    facet_wrap(~ Plot_ID, ncol = 3, nrow = 5, scales = "free_y") +
    labs(
      title = paste(title_prefix, "— Model_LAI (replicates 3×5)"),
      x = "Date", y = expression(~(mu*mol~m^{-2}~s^{-1}))
    ) +
    scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
          legend.position = "top")
}

# One-page PNG with all replicates (Calibration set)
if (nrow(cal_byrep)) {
  png("_Observed_vs_Modeled_AllReplicates_Model_LAI_BY_REP.png",
      width = 12, height = 15, units = "in", res = 300)
  print(timeseries_byrep(cal_byrep, "Observed vs Modeled – All 15 Replicates (Calibration)"))
  dev.off()
}

# Two-page PDFs (Calibration & Evaluation) — replicate panels on both pages
make_two_page_pdf_byrep <- function(df, file_out, title_prefix) {
  p1 <- timeseries_byrep(df, paste0(title_prefix, " (Page 1)"))
  p2 <- timeseries_byrep(df, paste0(title_prefix, " (Page 2)"))
  pdf(file_out, width = 14, height = 10)
  print(p1); print(p2)
  dev.off()
  message("✅ Saved: ", file_out)
}

if (nrow(cal_byrep))  make_two_page_pdf_byrep(cal_byrep,  "_Observed_vs_Modeled_CALIBRATION_Model_LAI_2pages_BY_REP.pdf", "Calibration")
if (nrow(eval_byrep)) make_two_page_pdf_byrep(eval_byrep, "_Observed_vs_Modeled_EVALUATION_Model_LAI_2pages_BY_REP.pdf",  "Evaluation")

# ---------- RMSE vs R² (BY REPLICATE) ----------
cal_metrics_rep <- cal_byrep %>%
  group_by(Veg, Rep) %>%
  summarise(
    RMSE_cal = hydroGOF::rmse(GPP_obs, GPP_mod),
    R2_cal   = { m <- lm(GPP_obs ~ GPP_mod); summary(m)$adj.r.squared },
    .groups = "drop"
  ) %>% mutate(Model = "Model_LAI")

eval_metrics_rep <- eval_byrep %>%
  group_by(Veg, Rep) %>%
  summarise(
    RMSE_eval = hydroGOF::rmse(GPP_obs, GPP_mod),
    R2_eval   = { m <- lm(GPP_obs ~ GPP_mod); summary(m)$adj.r.squared },
    .groups = "drop"
  ) %>% mutate(Model = "Model_LAI")

plot_cal_rep  <- ggplot(cal_metrics_rep,  aes(x = R2_cal,  y = RMSE_cal,  label = paste0(Veg, Rep))) +
  geom_point(color = "black", size = 3) +
  geom_text_repel(size = 3.2, max.overlaps = 100, box.padding = 0.3) +
  facet_wrap(~ Model) +
  scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 10), breaks = scales::pretty_breaks(n = 6)) +
  theme_bw(base_size = 14) +
  labs(title = "Calibration RMSE vs Adjusted R² — Model_LAI (by Replicate)",
       x = expression(Adjusted~R^2), y = "RMSE") +
  theme(strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5), panel.grid.major = element_line(color = "gray85"))

plot_eval_rep <- ggplot(eval_metrics_rep, aes(x = R2_eval, y = RMSE_eval, label = paste0(Veg, Rep))) +
  geom_point(color = "black", size = 3) +
  geom_text_repel(size = 3.2, max.overlaps = 100, box.padding = 0.3) +
  facet_wrap(~ Model) +
  scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 10), breaks = scales::pretty_breaks(n = 6)) +
  theme_bw(base_size = 14) +
  labs(title = "Evaluation RMSE vs Adjusted R² — Model_LAI (by Replicate)",
       x = expression(Adjusted~R^2), y = "RMSE") +
  theme(strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5), panel.grid.major = element_line(color = "gray85"))

combined_rmse_r2_rep <- (plot_cal_rep / plot_eval_rep) +
  patchwork::plot_annotation(title = "RMSE vs Adjusted R² — Model_LAI (Calibration top, Evaluation bottom; by Replicate)")

ggsave("RMSE_vs_R2_Calibration_and_Evaluation_Model_LAI_BY_REP.png",
       plot = combined_rmse_r2_rep, width = 12, height = 14, dpi = 300)

# ============================================================
#  REPLICATE-ONLY METRICS (Cal / Eval / Combined) — BY TREATMENT FIT
#  Outputs: GPP_Metrics_byRep_95CI_OE_UE_BY_TREATMENT.xlsx
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(openxlsx); library(stringr); library(hydroGOF)
})

# ---------- Helpers ----------
adj_r2_from_vectors <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) < 3) return(NA_real_)
  summary(lm(y[keep] ~ yhat[keep]))$adj.r.squared
}

# ---------- Build AIC map at Veg–Rep (fits are by Veg; we duplicate AIC to 1..3) ----------
if (exists("results_all_repl")) {
  aic_map <- results_all_repl %>%
    mutate(Rep = as.character(Rep), Model = "Model_LAI") %>%
    group_by(Veg, Rep, Model) %>%
    summarise(AIC = mean(AIC, na.rm = TRUE), .groups = "drop")
} else {
  aic_map <- results_all %>%
    mutate(Model = "Model_LAI") %>%
    tidyr::crossing(Rep = c("1","2","3")) %>%
    mutate(Rep = as.character(Rep)) %>%
    group_by(Veg, Rep, Model) %>%
    summarise(AIC = mean(AIC, na.rm = TRUE), .groups = "drop")
}

# ---------- Safety: ensure cal_byrep / eval_byrep exist ----------
if (!exists("cal_byrep"))  cal_byrep  <- data.frame()
if (!exists("eval_byrep")) eval_byrep <- data.frame()

# ---------- Metric calculator (per Veg–Rep–Model) ----------
compute_metrics_by_rep <- function(df, aic_map) {
  if (nrow(df) == 0) return(df %>% slice(0))
  df %>%
    group_by(Veg, Rep, Model) %>%
    summarise(
      n    = sum(is.finite(GPP_obs) & is.finite(GPP_mod)),
      RMSE = hydroGOF::rmse(GPP_obs, GPP_mod),
      R2_adj = adj_r2_from_vectors(GPP_obs, GPP_mod),
      OE_pct = mean(GPP_mod > GPP_obs, na.rm = TRUE) * 100,
      UE_pct = mean(GPP_mod < GPP_obs, na.rm = TRUE) * 100,
      Outside95_pct = {
        res <- GPP_obs - GPP_mod
        s <- sqrt(mean(res^2, na.rm = TRUE))          # group RMSE
        mean(abs(res) > 1.96 * s, na.rm = TRUE) * 100 # % outside ±1.96*RMSE
      },
      .groups = "drop"
    ) %>%
    left_join(aic_map, by = c("Veg","Rep","Model")) %>%
    arrange(factor(Veg, levels = c("B","BG","BGS","GS","G")), Rep)
}

# ---------- Compute per-replicate metrics ----------
cal_metrics_byrep  <- compute_metrics_by_rep(cal_byrep,  aic_map)
eval_metrics_byrep <- compute_metrics_by_rep(eval_byrep, aic_map)

# Combined = recompute metrics on all (Cal + Eval) points per Veg–Rep
combined_points <- bind_rows(cal_byrep, eval_byrep) %>% select(Veg, Rep, Model, GPP_obs, GPP_mod)
comb_metrics_byrep <- compute_metrics_by_rep(combined_points, aic_map)

# ---------- Save to Excel ----------
wb <- createWorkbook()
addWorksheet(wb, "Calibration_byRep"); writeData(wb, "Calibration_byRep", cal_metrics_byrep)
addWorksheet(wb, "Evaluation_byRep");  writeData(wb, "Evaluation_byRep",  eval_metrics_byrep)
addWorksheet(wb, "Combined_byRep");    writeData(wb, "Combined_byRep",    comb_metrics_byrep)
saveWorkbook(wb, "GPP_Metrics_byRep_95CI_OE_UE_BY_TREATMENT.xlsx", overwrite = TRUE)

cat("✅ Saved: GPP_Metrics_byRep_95CI_OE_UE_BY_TREATMENT.xlsx (replicate-only, Cal/Eval/Combined)\n")

