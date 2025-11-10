# =============================================================
# PART 1 — Base Setup and Reusable Functions
# =============================================================

# =============================
#      Load Required Libraries
# =============================
suppressPackageStartupMessages({
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
})

# =============================
#        Load and Prepare Data
# =============================
cbd <- read.csv("C:/Users/rboumelh/Desktop/Mesocosms/Modeling Reco/mesocosm_data1.csv", sep = ",")
cbd <- dplyr::rename(cbd, Date = TIMESTAMP.Reco)
cbd$Date <- as.POSIXct(cbd$Date, format = "%m/%d/%Y %H:%M")

# calibration / evaluation campaigns as before
CH4_cal  <- subset(cbd, ID_camp %in% c("1","3","4","6","8","10","11","13","15","16","18","19","21","22","24","25","27","28","30"))
CH4_eval <- subset(cbd, ID_camp %in% c("2","5","7","9","12","14","17","20","23","26"))

# we keep Tref (not used by the new model, but no harm to keep the variable)
Tref <- 288.15
setwd("C:/Users/rboumelh/Desktop/Mesocosms/Methane/Modeling CH4/3 temps_by treatment/other model/for git")

veg_types   <- c("B","BG","BGS","GS","G")
rep_levels  <- c("1","2","3")

# =============================
#   Generic Model Helper
#   NEW MODEL:
#   CH4 = (((a * (WTD / -0.35) + b)) + (c * LAI)) * ((T + 5)/20)^d
# =============================

predict_manual <- function(fit, newdata, temp_var, Tref = 288.15) {
  cf <- coef(fit)
  # make sure all vars exist
  Tval <- as.numeric(newdata[[temp_var]])
  WTDv <- as.numeric(newdata$WTD)
  LAIv <- as.numeric(newdata$LAI)
  
  # piece 1: linear eco part
  eco_part <- ((cf["a"] * (WTDv / -0.35)) + cf["b"]) + (cf["c"] * LAIv)
  
  # piece 2: temperature scaling
  temp_part <- ((Tval + 5) / 20) ^ cf["d"]
  
  eco_part * temp_part
}

needed_vars_for_model <- function(temp_var) {
  c("CH4", temp_var, "WTD", "LAI")
}

# =============================
#   Fitting Function (by vegetation)
#   → rewritten for the NEW MODEL
# =============================

fit_models_for_veg <- function(veg_name, CH4_cal, temp_var, model_label,
                               num_iter = 1, Tref = 288.15) {
  veg_data <- subset(CH4_cal, Plot == veg_name)
  out_results <- data.frame()
  out_preds   <- data.frame()
  fits        <- list()
  
  for (i in seq_len(num_iter)) {
    tryCatch({
      need_vars <- needed_vars_for_model(temp_var)
      rep_data_f <- veg_data %>%
        mutate(across(all_of(need_vars), as.numeric)) %>%
        filter(if_all(all_of(need_vars), ~ !is.na(.)))
      if (nrow(rep_data_f) == 0) stop("No valid rows for fitting")
      
      # NEW model formula
      formula_model <- as.formula(
        paste0(
          "CH4 ~ (((a * (WTD / -0.35) + b)) + (c * LAI)) * ((",
          temp_var,
          " + 5) / 20)^d"
        )
      )
      
      # starting values — conservative
      start_list   <- list(a = -1, b = 1, c = 1, d = 1)
      lower_bounds <- c(a = -Inf, b = -Inf, c = -Inf, d = -Inf)
      upper_bounds <- c(a =  Inf, b =  Inf, c =  Inf, d = Inf)
      
      fit <- nlsLM(
        formula_model, data = rep_data_f,
        start = start_list,
        lower = lower_bounds, upper = upper_bounds,
        control = nls.lm.control(maxiter = 1000, ftol = 1e-10)
      )
      
      # predictions on calibration rows
      preds    <- predict_manual(fit, rep_data_f, temp_var = temp_var, Tref = Tref)
      adj_r2   <- summary(lm(rep_data_f$CH4 ~ preds))$adj.r.squared
      rmse_val <- hydroGOF::rmse(rep_data_f$CH4, preds)
      aic_val  <- AIC(fit)
      bic_val  <- BIC(fit)
      
      model_key <- paste0(veg_name, "_iter", i, "_", model_label)
      fits[[model_key]] <- fit
      
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
        Model      = model_label,
        AIC        = aic_val,
        BIC        = bic_val,
        RMSE       = rmse_val,
        R2         = adj_r2,
        Veg        = veg_name
      )
      for (j in seq_along(params)) {
        pname <- rownames(coefs)[j]
        row[[pname]]                 <- params[j]
        row[[paste0("SD_", pname)]]  <- errors[j]
        row[[paste0("pval_", pname)]]<- p_vals[j]
      }
      out_results <- dplyr::bind_rows(out_results, as.data.frame(row))
      
      # store calibration predictions
      out_preds <- dplyr::bind_rows(out_preds, data.frame(
        Veg      = veg_name,
        Rep      = NA_character_,
        Model    = model_label,
        Date     = rep_data_f$Date,
        CH4_obs  = rep_data_f$CH4,
        CH4_mod  = preds,
        Residual = rep_data_f$CH4 - preds
      ))
    }, error = function(e) {
      cat("❌ Error in", veg_name, model_label, ":", e$message, "\n")
    })
  }
  
  list(results = out_results, preds = out_preds, fits = fits)
}

# =============================
#   Metrics Helper (unchanged)
# =============================
adj_r2_from_vectors <- function(y, yhat) {
  y    <- as.numeric(y)
  yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) < 3) return(NA_real_)
  summary(lm(y[keep] ~ yhat[keep]))$adj.r.squared
}

# =============================================================
# PART 2 — Full Dual-Model Workflow (Ts & Ta) with NEW MODEL
# =============================================================

R_gas <- 8.314   # left for compatibility
veg_levels_order <- c("B","BG","BGS","GS","G")
rep_levels_15    <- c("B1","B2","B3","BG1","BG2","BG3",
                      "BGS1","BGS2","BGS3",
                      "GS1","GS2","GS3",
                      "G1","G2","G3")

# -----------------------------
# One full end-to-end workflow
# -----------------------------
run_full_workflow <- function(temp_var, model_label) {
  
  message("\n\n========================")
  message(" Running ", model_label, " with ", temp_var)
  message("========================")
  
  # -----------------------------
  # 1) Fit by vegetation (plot-level)
  # -----------------------------
  all_res    <- list()
  all_preds  <- list()
  all_fits   <- list()
  
  for (veg in veg_types) {
    r <- fit_models_for_veg(veg, CH4_cal, temp_var, model_label, num_iter = 1, Tref = Tref)
    all_res[[veg]]   <- r$results
    all_preds[[veg]] <- r$preds
    all_fits         <- c(all_fits, r$fits)
  }
  
  results_all     <- dplyr::bind_rows(all_res)
  CH4_predictions <- dplyr::bind_rows(all_preds)
  fitted_models   <- all_fits
  
  # no Ea, no Q10 here — just order
  results_all <- results_all %>%
    mutate(
      Veg = factor(Veg, levels = veg_levels_order)
    )
  
  # Save parameters
  openxlsx::write.xlsx(results_all,
                       paste0("CH4_Model_Results_All_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # -----------------------------
  # 2) Calibration metrics (by replicate)
  # -----------------------------
  cal_results <- list()
  for (veg in veg_types) {
    for (rep_i in c("1","2","3")) {
      df_cal <- subset(CH4_cal, Plot == veg & Rep == rep_i)
      if (!nrow(df_cal)) next
      df_cal <- df_cal %>%
        mutate(across(all_of(needed_vars_for_model(temp_var)), as.numeric)) %>%
        filter(if_all(all_of(needed_vars_for_model(temp_var)), ~ is.finite(.)))
      if (!nrow(df_cal)) next
      
      fit <- fitted_models[[paste0(veg, "_iter1_", model_label)]]
      if (is.null(fit)) next
      
      preds <- predict_manual(fit, df_cal, temp_var = temp_var, Tref = Tref)
      keep  <- is.finite(df_cal$CH4) & is.finite(preds)
      if (sum(keep) < 3) next
      
      y_obs <- df_cal$CH4[keep]
      y_mod <- preds[keep]
      lin   <- lm(y_obs ~ y_mod)
      cal_results[[paste0(veg, rep_i)]] <- data.frame(
        Veg      = veg,
        Rep      = rep_i,
        Plot_ID  = paste0(veg, rep_i),
        Model    = model_label,
        RMSE_cal = hydroGOF::rmse(y_obs, y_mod),
        R2_cal   = summary(lin)$adj.r.squared,
        P_cal    = broom::glance(lin)$p.value
      )
    }
  }
  cal_results <- bind_rows(cal_results)
  openxlsx::write.xlsx(cal_results,
                       paste0("CH4_Model_Calibration_Metrics_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # -----------------------------
  # 3) Evaluation metrics (by replicate)
  # -----------------------------
  eval_results <- list()
  for (veg in veg_types) {
    for (rep_i in c("1","2","3")) {
      df_eval <- subset(CH4_eval, Plot == veg & Rep == rep_i)
      if (!nrow(df_eval)) next
      df_eval <- df_eval %>%
        mutate(across(all_of(needed_vars_for_model(temp_var)), as.numeric)) %>%
        filter(if_all(all_of(needed_vars_for_model(temp_var)), ~ is.finite(.)))
      if (!nrow(df_eval)) next
      
      fit <- fitted_models[[paste0(veg, "_iter1_", model_label)]]
      if (is.null(fit)) next
      
      preds <- predict_manual(fit, df_eval, temp_var = temp_var, Tref = Tref)
      keep  <- is.finite(df_eval$CH4) & is.finite(preds)
      if (sum(keep) < 3) next
      
      y_obs <- df_eval$CH4[keep]
      y_mod <- preds[keep]
      lin   <- lm(y_obs ~ y_mod)
      
      eval_results[[paste0(veg, rep_i)]] <- data.frame(
        Veg        = veg,
        Rep        = rep_i,
        Plot_ID    = paste0(veg, rep_i),
        Model      = model_label,
        RMSE_eval  = hydroGOF::rmse(y_obs, y_mod),
        R2_eval    = summary(lin)$adj.r.squared,
        P_eval     = broom::glance(lin)$p.value
      )
    }
  }
  eval_results <- bind_rows(eval_results)
  openxlsx::write.xlsx(eval_results,
                       paste0("CH4_Model_Evaluation_Metrics_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # -----------------------------
  # 4) Join perf (cal vs eval)
  # -----------------------------
  cal_tidy <- results_all %>%
    mutate(Rep = NA_character_) %>%
    select(Veg, Rep, Model,
           RMSE_cal = RMSE,
           R2_cal   = R2,
           AIC, BIC)
  
  perf_all <- eval_results %>%
    left_join(cal_tidy, by = c("Veg","Rep","Model"))
  openxlsx::write.xlsx(perf_all,
                       paste0("CH4_Model_Performance_Cal_vs_Eval_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # -----------------------------
  # 5) Scatter: Calibration/Evaluation by plot
  # -----------------------------
  CH4_predictions$Set <- "Calibration"
  
  # build eval predictions for saving
  eval_preds_df <- list()
  for (veg in veg_types) {
    fit <- fitted_models[[paste0(veg, "_iter1_", model_label)]]
    if (is.null(fit)) next
    subset_eval <- CH4_eval %>%
      filter(Plot == veg) %>%
      mutate(across(all_of(needed_vars_for_model(temp_var)), as.numeric)) %>%
      filter(if_all(all_of(needed_vars_for_model(temp_var)), ~ !is.na(.)))
    if (!nrow(subset_eval)) next
    preds_eval <- predict_manual(fit, subset_eval, temp_var = temp_var, Tref = Tref)
    eval_preds_df[[veg]] <- data.frame(
      Veg      = veg,
      Rep      = NA_character_,
      Model    = model_label,
      Date     = subset_eval$Date,
      CH4_obs  = subset_eval$CH4,
      CH4_mod  = preds_eval,
      Residual = subset_eval$CH4 - preds_eval,
      Set      = "Evaluation",
      stringsAsFactors = FALSE
    )
  }
  eval_preds_df <- bind_rows(eval_preds_df)
  
  make_sheet <- function(df) {
    if (nrow(df) == 0) return(df)
    df %>%
      mutate(Plot_ID = Veg,
             Date    = format(as.POSIXct(Date), "%m/%d/%Y %H:%M")) %>%
      select(Veg, Rep, Plot_ID, Model, Set, Date, CH4_obs, CH4_mod, Residual)
  }
  
  cal_to_save  <- make_sheet(CH4_predictions)
  eval_to_save <- make_sheet(eval_preds_df)
  comb_to_save <- make_sheet(bind_rows(CH4_predictions, eval_preds_df))
  
  openxlsx::write.xlsx(
    x = list(Calibration = cal_to_save,
             Evaluation  = eval_to_save,
             Combined    = comb_to_save),
    file = paste0("CH4_Observed_vs_Modeled_All_", model_label, ".xlsx"),
    overwrite = TRUE
  )
  
  # Facet order
  combined_preds <- bind_rows(CH4_predictions, eval_preds_df) %>%
    mutate(Plot_ID = factor(Veg, levels = veg_levels_order))
  
  safe_stats <- function(df) {
    ok <- is.finite(df$CH4_obs) & is.finite(df$CH4_mod)
    df <- df[ok, , drop = TRUE]
    if (nrow(df) < 2 || stats::sd(df$CH4_obs) == 0) {
      return(tibble(R2 = NA_real_, Pval = NA_real_))
    }
    m <- lm(CH4_mod ~ CH4_obs, data = df)
    tibble(R2 = summary(m)$adj.r.squared,
           Pval = broom::glance(m)$p.value)
  }
  
  stats_df <- combined_preds %>%
    group_split(Plot_ID, Set, Model) %>%
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
    }) %>% bind_rows()
  
  openxlsx::write.xlsx(stats_df,
                       paste0("CH4_Panel_Stats_R2_Pval_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  combined_preds <- left_join(combined_preds, stats_df,
                              by = c("Plot_ID","Set","Model"))
  
  plot_scatter_by_set <- function(set_type) {
    df <- combined_preds %>% filter(Set == set_type, Model == model_label)
    annotation_df <- df %>%
      group_by(Plot_ID) %>%
      slice(1) %>%
      mutate(label = paste0("R² = ", round(R2, 2), "\nP = ", signif(Pval, 2)))
    ggplot(df, aes(x = CH4_obs, y = CH4_mod)) +
      geom_point(alpha = 0.6, color = "black", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, color = "gray40", linetype = "dashed") +
      geom_text(data = annotation_df,
                aes(x = Inf, y = Inf, label = label),
                inherit.aes = FALSE, size = 4.5, hjust = 1.1, vjust = 1.5) +
      facet_wrap(~ Plot_ID, ncol = 3, scales = "free") +
      labs(
        title = paste(set_type, "-", model_label, "- Observed vs Modeled CH4 (by Plot)"),
        x = expression(Observed~CH4~(mu*mol~m^{-2}~s^{-1})),
        y = expression(Modeled~CH4~(mu*mol~m^{-2}~s^{-1}))
      ) +
      theme_bw(base_size = 12) +
      theme(strip.text = element_text(size = 11, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            axis.text = element_text(size = 9),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank())
  }
  
  fig_cal  <- plot_scatter_by_set("Calibration")
  fig_eval <- plot_scatter_by_set("Evaluation")
  ggsave(paste0("CH4_Calibration_Scatterplots_byPlot_", model_label, ".png"),
         fig_cal,  width = 10, height = 7, dpi = 300)
  ggsave(paste0("CH4_Evaluation_Scatterplots_byPlot_",  model_label, ".png"),
         fig_eval, width = 10, height = 7, dpi = 300)
  
  # -----------------------------
  # 6) Time series (calibration) 3×2 — keep single PNG only
  # -----------------------------
  data_ts <- CH4_predictions %>%
    mutate(Date = as.POSIXct(Date),
           Plot_ID = factor(Veg, levels = veg_levels_order)) %>%
    pivot_longer(cols = c(CH4_obs, CH4_mod),
                 names_to = "Type", values_to = "CH4") %>%
    filter(!is.na(CH4)) %>%
    mutate(Type = recode(Type,
                         CH4_obs = "Observed",
                         CH4_mod = "Modeled"))
  
  n_panels  <- length(unique(data_ts$Plot_ID))
  ncol_auto <- 3
  nrow_auto <- ceiling(n_panels / ncol_auto)
  
  g_ts_all <- ggplot(data_ts, aes(x = Date, y = CH4, color = Type)) +
    geom_point(size = 1.4) +
    facet_wrap(~ Plot_ID, ncol = ncol_auto, nrow = nrow_auto, scales = "free_y") +
    labs(title = paste("Observed vs Modeled – All Plots (", model_label, ")", sep = ""),
         x = "Date",
         y = expression(~(mu*mol~m^{-2}~s^{-1}))) +
    scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
    theme_bw(base_size = 14) +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 90, size = 9),
          axis.text.y = element_text(size = 9),
          legend.position = "top")
  ggsave(paste0("_Observed_vs_Modeled_AllPlots_", model_label, ".png"),
         g_ts_all, width = 12, height = 8, dpi = 300)
  
  # (Paginated PDFs removed per request)
  
  # -----------------------------
  # 7) Build replicate-level predictions (B1..G3) from plot-level fit
  # -----------------------------
  .parse_key_to_veg <- function(model_key) {
    mi <- strcapture("^([A-Z]+)_iter[0-9]+_(.+)$", model_key,
                     proto = list(Veg = character(), Model = character()))
    mi$Veg
  }
  
  build_replicate_preds <- function(base_df, set_name) {
    out <- list()
    for (model_key in names(fitted_models)) {
      fit  <- fitted_models[[model_key]]
      veg  <- .parse_key_to_veg(model_key)
      need <- needed_vars_for_model(temp_var)
      for (rep_i in c("1","2","3")) {
        df <- base_df %>%
          filter(Plot == veg, Rep == rep_i) %>%
          mutate(across(all_of(need), as.numeric)) %>%
          filter(if_all(all_of(need), ~ is.finite(.)))
        if (!nrow(df)) next
        preds <- predict_manual(fit, df, temp_var = temp_var, Tref = Tref)
        out[[length(out) + 1]] <- data.frame(
          Veg      = veg,
          Rep      = rep_i,
          Plot_ID  = paste0(veg, rep_i),
          Model    = model_label,
          Date     = df$Date,
          CH4_obs  = df$CH4,
          CH4_mod  = preds,
          Residual = df$CH4 - preds,
          Set      = set_name,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(out)) bind_rows(out) else
      data.frame(Veg = character(), Rep = character(), Plot_ID = character(),
                 Model = character(), Date = as.POSIXct(character()),
                 CH4_obs = numeric(), CH4_mod = numeric(), Residual = numeric(),
                 Set = character(), stringsAsFactors = FALSE)
  }
  
  cal_byrep  <- build_replicate_preds(CH4_cal,  "Calibration")
  eval_byrep <- build_replicate_preds(CH4_eval, "Evaluation")
  
  # -----------------------------
  # 8) SCATTERPLOTS (BY REPLICATE)
  # -----------------------------
  scatter_byrep <- function(df, title_txt) {
    stats_df <- df %>%
      group_by(Plot_ID) %>%
      summarise(
        R2 = {
          m <- try(lm(CH4_mod ~ CH4_obs, data = .), silent = TRUE)
          if (inherits(m, "try-error")) NA_real_ else summary(m)$adj.r.squared
        },
        Pval = {
          m <- try(lm(CH4_mod ~ CH4_obs, data = .), silent = TRUE)
          if (inherits(m, "try-error")) NA_real_ else broom::glance(m)$p.value
        },
        .groups = "drop"
      )
    df2 <- df %>% mutate(Plot_ID = factor(Plot_ID, levels = rep_levels_15))
    ann <- stats_df %>%
      mutate(label = paste0("R² = ", round(R2, 2), "\nP = ", signif(Pval, 2)))
    ggplot(df2, aes(x = CH4_obs, y = CH4_mod)) +
      geom_point(alpha = 0.6, size = 1.2, color = "black") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
      geom_text(data = ann, aes(x = Inf, y = Inf, label = label),
                inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 4.2) +
      facet_wrap(~ Plot_ID, ncol = 3, scales = "free") +
      labs(title = title_txt,
           x = expression(Observed~CH4~(mu*mol~m^{-2}~s^{-1})),
           y = expression(Modeled~CH4~(mu*mol~m^{-2}~s^{-1}))) +
      theme_bw(base_size = 12) +
      theme(strip.text = element_text(size = 11, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank())
  }
  
  if (nrow(cal_byrep)) {
    p_cal_rep <- scatter_byrep(
      cal_byrep,
      paste("Calibration —", model_label, "— Observed vs Modeled (by Replicate)")
    )
    ggsave(paste0("CH4_Calibration_Scatterplots_byReplicate_", model_label, ".png"),
           p_cal_rep, width = 10, height = 15, dpi = 300)
  }
  if (nrow(eval_byrep)) {
    p_eval_rep <- scatter_byrep(
      eval_byrep,
      paste("Evaluation —", model_label, "— Observed vs Modeled (by Replicate)")
    )
    ggsave(paste0("CH4_Evaluation_Scatterplots_byReplicate_", model_label, ".png"),
           p_eval_rep, width = 10, height = 15, dpi = 300)
  }
  
  # -----------------------------
  # 9) TIME SERIES 3×5 (BY REPLICATE)
  # -----------------------------
  timeseries_byrep <- function(df, title_prefix) {
    df2 <- df %>%
      mutate(Plot_ID = factor(Plot_ID, levels = rep_levels_15)) %>%
      pivot_longer(c(CH4_obs, CH4_mod), names_to = "Type", values_to = "CH4") %>%
      mutate(Type = recode(Type, CH4_obs = "Observed", CH4_mod = "Modeled"))
    ggplot(df2, aes(x = Date, y = CH4, color = Type)) +
      geom_point(size = 1.2, na.rm = TRUE) +
      facet_wrap(~ Plot_ID, ncol = 3, scales = "free_y") +
      labs(
        title = paste(title_prefix, "—", model_label, "(replicates 3×5)"),
        x = "Date", y = expression(~(mu*mol~m^{-2}~s^{-1}))
      ) +
      scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
      theme_bw(base_size = 13) +
      theme(strip.text = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
            legend.position = "top")
  }
  
  if (nrow(cal_byrep)) {
    png(paste0("_Observed_vs_Modeled_AllReplicates_", model_label, ".png"),
        width = 12, height = 15, units = "in", res = 300)
    print(timeseries_byrep(cal_byrep, "Observed vs Modeled – All 15 Replicates (Calibration)"))
    dev.off()
  }
  
  # -----------------------------
  # 10) RMSE vs R² — Plot-level
  # -----------------------------
  plot_cal <- ggplot(cal_results, aes(x = R2_cal, y = RMSE_cal, label = Veg)) +
    geom_point(color = "black", size = 3) +
    geom_text_repel(size = 3.2, max.overlaps = 100, box.padding = 0.3) +
    facet_wrap(~ Model, scales = "fixed") +
    scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(breaks = pretty_breaks(n = 6)) +
    theme_bw(base_size = 14) +
    labs(title = paste("Calibration RMSE vs Adjusted R² —", model_label),
         x = expression(Adjusted~R^2),
         y = "RMSE") +
    theme(strip.text = element_text(size = 14, face = "bold"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_line(color = "gray85"))
  ggsave(paste0("Calibration_RMSE_vs_R2_", model_label, ".png"),
         plot_cal, width = 12, height = 8, dpi = 300)
  
  plot_eval <- ggplot(eval_results,
                      aes(x = R2_eval, y = RMSE_eval, label = Veg)) +
    geom_point(color = "black", size = 3) +
    ggrepel::geom_text_repel(size = 3.2, max.overlaps = 100, box.padding = 0.3) +
    facet_wrap(~ Model, scales = "fixed") +
    scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    theme_bw(base_size = 14) +
    labs(title = paste("Evaluation RMSE vs Adjusted R² —", model_label),
         x = expression(Adjusted~R^2),
         y = "RMSE") +
    theme(strip.text = element_text(size = 14, face = "bold"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_line(color = "gray85"))
  ggsave(paste0("Evaluation_RMSE_vs_R2_", model_label, ".png"),
         plot = plot_eval, width = 12, height = 8, dpi = 300)
  
  # -----------------------------
  # 11) Parameter histograms & summaries (a, b, c, d)
  # -----------------------------
  params <- results_all %>%
    mutate(
      SubsetType = Veg
    )
  
  subset_levels <- veg_levels_order
  
  params_tidy <- params %>%
    transmute(
      SubsetType, Model,
      a, b, c, d,
      SD_a = dplyr::coalesce(SD_a, NA_real_),
      SD_b = dplyr::coalesce(SD_b, NA_real_),
      SD_c = dplyr::coalesce(SD_c, NA_real_),
      SD_d = dplyr::coalesce(SD_d, NA_real_),
      pval_a = dplyr::coalesce(pval_a, NA_real_),
      pval_b = dplyr::coalesce(pval_b, NA_real_),
      pval_c = dplyr::coalesce(pval_c, NA_real_),
      pval_d = dplyr::coalesce(pval_d, NA_real_)
    ) %>%
    pivot_longer(cols = c(a, b, c, d),
                 names_to = "Parameter",
                 values_to = "Estimate") %>%
    mutate(
      SE = case_when(
        Parameter == "a" ~ SD_a,
        Parameter == "b" ~ SD_b,
        Parameter == "c" ~ SD_c,
        Parameter == "d" ~ SD_d
      ),
      Pval = case_when(
        Parameter == "a" ~ pval_a,
        Parameter == "b" ~ pval_b,
        Parameter == "c" ~ pval_c,
        Parameter == "d" ~ pval_d
      ),
      Signif     = ifelse(!is.na(Pval) & Pval < 0.05, "*", ""),
      SubsetType = factor(SubsetType, levels = subset_levels),
      Model      = factor(Model)
    )
  
  plot_param_hist <- function(param_name, ylab_default = "Estimate ± SE") {
    df <- filter(params_tidy, Parameter == param_name)
    if (nrow(df) == 0) {
      warning(paste("No data for", param_name))
      return(NULL)
    }
    ggplot(df, aes(x = SubsetType, y = Estimate)) +
      geom_col(color = "black", fill = NA, width = 0.7) +
      geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2) +
      geom_text(aes(label = Signif, y = Estimate + SE + 0.1),
                vjust = 0, size = 5, na.rm = TRUE) +
      facet_grid(cols = vars(Model), scales = "free_y") +
      labs(title = paste("Parameter ", param_name, sep = ""),
           y = ylab_default, x = "Plot") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", size = 16),
            strip.text = element_text(size = 12, face = "bold"),
            legend.position = "none")
  }
  
  pdf(paste0("Parameter_a_b_c_d_Histograms_", model_label, ".pdf"),
      width = 14, height = 8)
  for (p in c("a","b","c","d")) {
    p_obj <- plot_param_hist(p)
    if (!is.null(p_obj)) print(p_obj)
  }
  dev.off()
  
  # Group summaries (means ± SE)
  long_params <- params %>%
    transmute(
      Group     = Veg,
      Replicate = Veg,
      Model,
      a, b, c, d,
      SE_a = dplyr::coalesce(SD_a, NA_real_),
      SE_b = dplyr::coalesce(SD_b, NA_real_),
      SE_c = dplyr::coalesce(SD_c, NA_real_),
      SE_d = dplyr::coalesce(SD_d, NA_real_),
      P_a  = dplyr::coalesce(pval_a, NA_real_),
      P_b  = dplyr::coalesce(pval_b, NA_real_),
      P_c  = dplyr::coalesce(pval_c, NA_real_),
      P_d  = dplyr::coalesce(pval_d, NA_real_)
    ) %>%
    pivot_longer(cols = c(a, b, c, d),
                 names_to = "Parameter",
                 values_to = "Estimate") %>%
    mutate(
      SE = case_when(
        Parameter == "a" ~ SE_a,
        Parameter == "b" ~ SE_b,
        Parameter == "c" ~ SE_c,
        Parameter == "d" ~ SE_d
      ),
      Pval = case_when(
        Parameter == "a" ~ P_a,
        Parameter == "b" ~ P_b,
        Parameter == "c" ~ P_c,
        Parameter == "d" ~ P_d
      )
    ) %>%
    select(Group, Replicate, Model, Parameter, Estimate, SE, Pval)
  
  summary_params <- long_params %>%
    group_by(Group, Model, Parameter) %>%
    summarise(
      Mean_Estimate = mean(Estimate, na.rm = TRUE),
      SE_Estimate   = sd(Estimate,   na.rm = TRUE) / sqrt(sum(!is.na(Estimate))),
      Min_Pval      = suppressWarnings(min(Pval, na.rm = TRUE)),
      .groups       = "drop"
    ) %>%
    mutate(
      Signif = case_when(
        is.finite(Min_Pval) & Min_Pval < 0.001 ~ "***",
        is.finite(Min_Pval) & Min_Pval < 0.01  ~ "**",
        is.finite(Min_Pval) & Min_Pval < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Model = factor(Model),
      Group = factor(Group, levels = veg_levels_order)
    )
  
  # (Averaged_Parameter_Histograms_* removed per request — no PDF saved)
  
  # -----------------------------
  # 12) Keep results in .GlobalEnv to feed other figures
  # -----------------------------
  assign(paste0("RESULTS_", model_label), results_all, envir = .GlobalEnv)
  assign(paste0("CH4CAL_",  model_label), CH4_cal,    envir = .GlobalEnv)
  assign(paste0("TEMPVAR_", model_label), temp_var,   envir = .GlobalEnv)
  
  message("✅ Finished: ", model_label)
  invisible(list(
    results_all     = results_all,
    CH4_predictions = CH4_predictions,
    eval_results    = eval_results,
    cal_results     = cal_results,
    cal_byrep       = cal_byrep,
    eval_byrep      = eval_byrep,
    fitted_models   = fitted_models
  ))
}

# -----------------------------
# Run the three models (Ts15, Ts5, Ts30) with NEW formula
# -----------------------------
out_Ts15 <- run_full_workflow("Ts15",    "Model_Ts15")
out_Ts5  <- run_full_workflow("RE_5_Ts", "Model_Ts5")
# Ts30 will be run later (as in your original code)

# =============================================================
# PART 3 — Facet Order Utilities (Force B, BG, BGS, GS, G)
# =============================================================
force_veg_order <- function(df, column = "Group") {
  veg_levels <- c("B","BG","BGS","GS","G")
  if (column %in% names(df)) {
    df[[column]] <- factor(df[[column]], levels = veg_levels)
  }
  df
}

reorder_facets_B_BG_BGS_GS_G <- function(p) {
  veg_levels <- c("B","BG","BGS","GS","G")
  p + facet_wrap(
    vars(Group),
    ncol = length(veg_levels),
    scales = "free_x",
    labeller = label_value
  ) +
    scale_x_discrete(limits = veg_levels)
}

# =============================================================
# COMBINED 2×2 PANEL — Calibration & Evaluation (Ts5 vs Ts15)
# =============================================================
library(patchwork)
library(ggplot2)
library(png)
library(grid)

plot_cal_Ts5  <- ggplot2::ggplot() +
  ggplot2::annotation_custom(
    grid::rasterGrob(png::readPNG("Calibration_RMSE_vs_R2_Model_Ts5.png")),
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  ggplot2::theme_void()

plot_eval_Ts5 <- ggplot2::ggplot() +
  ggplot2::annotation_custom(
    grid::rasterGrob(png::readPNG("Evaluation_RMSE_vs_R2_Model_Ts5.png")),
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  ggplot2::theme_void()

plot_cal_Ts15  <- ggplot2::ggplot() +
  ggplot2::annotation_custom(
    grid::rasterGrob(png::readPNG("Calibration_RMSE_vs_R2_Model_Ts15.png")),
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  ggplot2::theme_void()

plot_eval_Ts15 <- ggplot2::ggplot() +
  ggplot2::annotation_custom(
    grid::rasterGrob(png::readPNG("Evaluation_RMSE_vs_R2_Model_Ts15.png")),
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  ggplot2::theme_void()

combined_panel <- (plot_cal_Ts5 | plot_cal_Ts15) /
  (plot_eval_Ts5 | plot_eval_Ts15) +
  patchwork::plot_annotation(
    title = "Calibration and Evaluation RMSE vs Adjusted R² (NEW model)",
    subtitle = "Left: Ts5 model | Right: Ts15 model",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

ggsave("RMSE_vs_R2_Ts5_vs_Ts15_2x2.png",
       combined_panel, width = 14, height = 10, dpi = 300)
ggsave("RMSE_vs_R2_Ts5_vs_Ts15_2x2.pdf",
       combined_panel, width = 14, height = 10)

cat("✅ Saved: RMSE_vs_R2_Ts5_vs_Ts15_2x2.(png|pdf)\n")

# =============================================================
# 2×2 PANEL — RMSE vs R² BY REPLICATE (B1…G3): Ts5 vs Ts15
# =============================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel);
  library(scales); library(patchwork); library(hydroGOF)
})

.compute_rep_metrics <- function(df, set_label = c("Calibration","Evaluation")) {
  set_label <- match.arg(set_label)
  if (!nrow(df)) return(data.frame())
  df %>%
    group_by(Veg, Rep) %>%
    summarise(
      RMSE = hydroGOF::rmse(CH4_obs, CH4_mod),
      R2   = { m <- lm(CH4_obs ~ CH4_mod); summary(m)$adj.r.squared },
      .groups = "drop"
    ) %>%
    mutate(PanelLabel = paste0(Veg, Rep), Set = set_label)
}

# compute metrics
cal_rep_Ts5   <- .compute_rep_metrics(out_Ts5$cal_byrep,  "Calibration")
eval_rep_Ts5  <- .compute_rep_metrics(out_Ts5$eval_byrep, "Evaluation")
cal_rep_Ts15  <- .compute_rep_metrics(out_Ts15$cal_byrep,  "Calibration")
eval_rep_Ts15 <- .compute_rep_metrics(out_Ts15$eval_byrep, "Evaluation")

.plot_rmse_r2_byrep <- function(df, title_txt) {
  if (!nrow(df)) {
    return(ggplot() + ggtitle(paste(title_txt, "(no data)")) + theme_void())
  }
  
  df$Veg <- factor(df$Veg, levels = c("B","BG","BGS","GS","G"))
  df <- df %>% arrange(Veg, Rep)
  
  # no p-value coloring in this version
  df$SigColor <- "black"
  
  ggplot(df, aes(x = R2, y = RMSE, label = PanelLabel)) +
    geom_point(aes(color = SigColor), size = 2.8) +
    scale_color_identity() +
    geom_text_repel(size = 3.1, max.overlaps = 100,
                    box.padding = 0.3, min.segment.length = 0,
                    color = "black") +
    scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02)) +
    labs(title = title_txt,
         x = expression(Adjusted~R^2),
         y = "RMSE") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank())
}

# ============================
# RUN THIRD MODEL (Ts30) — NEW MODEL
# ============================
out_Ts30 <- run_full_workflow("Ts30", "Model_Ts30")

# build metrics for Ts30
cal_rep_Ts30  <- .compute_rep_metrics(out_Ts30$cal_byrep,  "Calibration")
eval_rep_Ts30 <- .compute_rep_metrics(out_Ts30$eval_byrep, "Evaluation")

# plots
p_cal_Ts5   <- .plot_rmse_r2_byrep(cal_rep_Ts5,   "Calibration — Ts5")
p_cal_Ts15  <- .plot_rmse_r2_byrep(cal_rep_Ts15,  "Calibration — Ts15")
p_cal_Ts30  <- .plot_rmse_r2_byrep(cal_rep_Ts30,  "Calibration — Ts30")

p_eval_Ts5  <- .plot_rmse_r2_byrep(eval_rep_Ts5,  "Evaluation — Ts5")
p_eval_Ts15 <- .plot_rmse_r2_byrep(eval_rep_Ts15, "Evaluation — Ts15")
p_eval_Ts30 <- .plot_rmse_r2_byrep(eval_rep_Ts30, "Evaluation — Ts30")

rmse_r2_byrep_3x2 <- (p_cal_Ts5 | p_cal_Ts15 | p_cal_Ts30) /
  (p_eval_Ts5 | p_eval_Ts15 | p_eval_Ts30) +
  patchwork::plot_annotation(
    title = "RMSE vs Adjusted R² by Replicate — Ts5 vs Ts15 vs Ts30 (NEW model)",
    subtitle = "Top: Calibration | Bottom: Evaluation",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

ggsave("RMSE_vs_R2_byRep_Ts5_Ts15_Ts30_3x2.png",
       rmse_r2_byrep_3x2, width = 20, height = 10, dpi = 300)
ggsave("RMSE_vs_R2_byRep_Ts5_Ts15_Ts30_3x2.pdf",
       rmse_r2_byrep_3x2, width = 20, height = 10)

cat("✅ Saved: RMSE_vs_R2_byRep_Ts5_Ts15_Ts30_3x2.(png|pdf)\n")
