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

# same splits as before
Reco_cal  <- subset(cbd, ID_camp %in% c("1","3","4","6","8","10","11","13","15","16","18","19","21","22","24","25","27","28","30"))
Reco_eval <- subset(cbd, ID_camp %in% c("2","5","7","9","12","14","17","20","23","26"))

setwd("C:/Users/rboumelh/Desktop/Mesocosms/Modeling Reco/Modeling/Ta&Ts_by treat/other model/for git")

veg_types  <- c("B","BG","BGS","GS","G")
rep_levels <- c("1","2","3")

# =============================
#   Generic Model Helpers
# =============================

# We now predict with the NEW model:
# Reco = (((a * (WTD / -0.35) + b)) + (c * LAI)) * ((T + 5) / 20)^d
# temp_var is either "RE_Ta" or "RE_5_Ts"
predict_manual <- function(fit, newdata, temp_var) {
  cf <- coef(fit)
  Tval <- newdata[[temp_var]]
  base_lin <- (cf["a"] * (newdata$WTD / -0.35)) + cf["b"] + (cf["c"] * newdata$LAI)
  temp_fac <- ((Tval + 5) / 20)^cf["d"]
  base_lin * temp_fac
}

needed_vars_for_model <- function(temp_var) {
  c("Reco", temp_var, "WTD", "LAI")
}

# =============================
#   Fitting Function (by vegetation)
#   → uses the NEW formula
# =============================
fit_models_for_veg <- function(veg_name, Reco_cal, temp_var, model_label, num_iter = 1) {
  veg_data <- subset(Reco_cal, Plot == veg_name)
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
      
      # ---- NEW MODEL ----
      # Reco ~ (((a * (WTD / -0.35) + b)) + (c * LAI)) * ((T + 5) / 20)^d
      formula_model <- as.formula(
        paste0(
          "Reco ~ (((a * (WTD / -0.35)) + b) + (c * LAI)) * ((",
          temp_var,
          " + 5) / 20)^d"
        )
      )
      
      # starting values — gentle
      start_list   <- list(a = -0.5, b = 0.01, c = 0.01, d = 1)
      lower_bounds <- c(a = -Inf, b = -Inf, c = -Inf, d = -5)
      upper_bounds <- c(a =  Inf, b =  Inf, c =  Inf, d =  5)
      
      fit <- nlsLM(
        formula_model,
        data   = rep_data_f,
        start  = start_list,
        lower  = lower_bounds,
        upper  = upper_bounds,
        control = nls.lm.control(maxiter = 1000, ftol = 1e-10)
      )
      
      preds    <- predict(fit)
      adj_r2   <- summary(lm(rep_data_f$Reco ~ preds))$adj.r.squared
      rmse_val <- rmse(rep_data_f$Reco, preds)
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
        row[[pname]]                   <- params[j]
        row[[paste0("SD_", pname)]]    <- errors[j]
        row[[paste0("pval_", pname)]]  <- p_vals[j]
      }
      out_results <- bind_rows(out_results, as.data.frame(row))
      
      out_preds <- bind_rows(out_preds, data.frame(
        Veg      = veg_name,
        Rep      = NA_character_,
        Model    = model_label,
        Date     = rep_data_f$Date,
        Reco_obs = rep_data_f$Reco,
        Reco_mod = preds,
        Residual = rep_data_f$Reco - preds
      ))
    }, error = function(e) {
      cat("❌ Error in", veg_name, model_label, ":", e$message, "\n")
    })
  }
  
  list(results = out_results, preds = out_preds, fits = fits)
}

# =============================
#   Metrics Helper
# =============================
adj_r2_from_vectors <- function(y, yhat) {
  y <- as.numeric(y); yhat <- as.numeric(yhat)
  keep <- is.finite(y) & is.finite(yhat)
  if (sum(keep) < 3) return(NA_real_)
  summary(lm(y[keep] ~ yhat[keep]))$adj.r.squared
}

# =============================================================
# PART 2 — Full Dual-Model Workflow (Ta & Ts)
# =============================================================

veg_levels_order <- c("B","BG","BGS","GS","G")
rep_levels_15    <- c("B1","B2","B3",
                      "BG1","BG2","BG3",
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
  
  # 1) Fit by vegetation
  all_res   <- list(); all_preds <- list(); all_fits <- list()
  for (veg in veg_types) {
    r <- fit_models_for_veg(veg, Reco_cal, temp_var, model_label, num_iter = 1)
    all_res[[veg]]   <- r$results
    all_preds[[veg]] <- r$preds
    all_fits         <- c(all_fits, r$fits)
  }
  results_all      <- bind_rows(all_res)
  Reco_predictions <- bind_rows(all_preds)
  fitted_models    <- all_fits
  
  # (no Ea, no Q10 now)
  results_all <- results_all %>%
    mutate(Veg = factor(Veg, levels = veg_levels_order))
  
  # Save parameters
  openxlsx::write.xlsx(results_all,
                       paste0("Reco_Model_Results_All_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # 2) Calibration metrics (plot-level)
  cal_results <- Reco_predictions %>%
    mutate(Veg = factor(Veg, levels = veg_levels_order)) %>%
    group_by(Veg, Model) %>%
    summarise(
      {
        keep <- is.finite(Reco_obs) & is.finite(Reco_mod)
        if (sum(keep) < 3) {
          tibble(RMSE_cal = NA_real_, R2_cal = NA_real_, P_cal = NA_real_)
        } else {
          y_obs <- Reco_obs[keep]; y_mod <- Reco_mod[keep]
          lin   <- lm(y_obs ~ y_mod)
          tibble(
            RMSE_cal = hydroGOF::rmse(y_obs, y_mod),
            R2_cal   = summary(lin)$adj.r.squared,
            P_cal    = broom::glance(lin)$p.value
          )
        }
      },
      .groups = "drop"
    ) %>%
    mutate(Rep = NA_character_) %>%
    relocate(Rep, .after = Veg)
  openxlsx::write.xlsx(cal_results,
                       paste0("Reco_Model_Calibration_Metrics_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # 3) Evaluation (plot-level)
  eval_results <- list()
  for (veg in veg_types) {
    df_eval <- subset(Reco_eval, Plot == veg) %>%
      mutate(across(all_of(needed_vars_for_model(temp_var)), as.numeric)) %>%
      filter(if_all(all_of(needed_vars_for_model(temp_var)), ~ is.finite(.)))
    if (!nrow(df_eval)) next
    fit <- fitted_models[[paste0(veg, "_iter1_", model_label)]]
    if (is.null(fit)) next
    preds_eval <- predict_manual(fit, df_eval, temp_var = temp_var)
    keep <- is.finite(df_eval$Reco) & is.finite(preds_eval)
    if (sum(keep) < 3) next
    y_obs <- df_eval$Reco[keep]; y_mod <- preds_eval[keep]
    lin   <- lm(y_obs ~ y_mod)
    eval_results[[veg]] <- data.frame(
      Veg       = veg,
      Rep       = NA_character_,
      Model     = model_label,
      R2_eval   = summary(lin)$adj.r.squared,
      RMSE_eval = hydroGOF::rmse(y_obs, y_mod),
      P_eval    = broom::glance(lin)$p.value
    )
  }
  eval_results <- bind_rows(eval_results)
  openxlsx::write.xlsx(eval_results,
                       paste0("Reco_Model_Evaluation_Metrics_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # 4) Join perf (cal vs eval)
  cal_tidy <- results_all %>%
    mutate(Rep = NA_character_) %>%
    select(Veg, Rep, Model,
           RMSE_cal = RMSE,
           R2_cal   = R2,
           AIC, BIC)
  perf_all <- eval_results %>%
    left_join(cal_tidy, by = c("Veg","Rep","Model"))
  openxlsx::write.xlsx(perf_all,
                       paste0("Reco_Model_Performance_Cal_vs_Eval_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  # 5) Save Obs/Mod (cal+eval) and scatterplots
  Reco_predictions$Set <- "Calibration"
  
  # build eval predictions for saving
  eval_preds_df <- list()
  for (veg in veg_types) {
    fit <- fitted_models[[paste0(veg, "_iter1_", model_label)]]
    if (is.null(fit)) next
    subset_eval <- Reco_eval %>%
      filter(Plot == veg) %>%
      mutate(across(all_of(needed_vars_for_model(temp_var)), as.numeric)) %>%
      filter(if_all(all_of(needed_vars_for_model(temp_var)), ~ !is.na(.)))
    if (!nrow(subset_eval)) next
    preds_eval <- predict_manual(fit, subset_eval, temp_var = temp_var)
    eval_preds_df[[veg]] <- data.frame(
      Veg      = veg,
      Rep      = NA_character_,
      Model    = model_label,
      Date     = subset_eval$Date,
      Reco_obs = subset_eval$Reco,
      Reco_mod = preds_eval,
      Residual = subset_eval$Reco - preds_eval,
      Set      = "Evaluation",
      stringsAsFactors = FALSE
    )
  }
  eval_preds_df <- bind_rows(eval_preds_df)
  
  make_sheet <- function(df) {
    if (nrow(df) == 0) return(df)
    df %>%
      mutate(Plot_ID = Veg,
             Date = format(as.POSIXct(Date), "%m/%d/%Y %H:%M")) %>%
      select(Veg, Rep, Plot_ID, Model, Set, Date, Reco_obs, Reco_mod, Residual)
  }
  cal_to_save  <- make_sheet(Reco_predictions)
  eval_to_save <- make_sheet(eval_preds_df)
  comb_to_save <- make_sheet(bind_rows(Reco_predictions, eval_preds_df))
  openxlsx::write.xlsx(
    x = list(Calibration = cal_to_save,
             Evaluation = eval_to_save,
             Combined   = comb_to_save),
    file = paste0("Reco_Observed_vs_Modeled_All_", model_label, ".xlsx"),
    overwrite = TRUE
  )
  
  # scatter by set
  combined_preds <- bind_rows(Reco_predictions, eval_preds_df) %>%
    mutate(Plot_ID = factor(Veg, levels = veg_levels_order))
  
  safe_stats <- function(df) {
    ok <- is.finite(df$Reco_obs) & is.finite(df$Reco_mod)
    df <- df[ok, , drop = TRUE]
    if (nrow(df) < 2 || stats::sd(df$Reco_obs) == 0) {
      return(tibble(R2 = NA_real_, Pval = NA_real_))
    }
    m <- lm(Reco_mod ~ Reco_obs, data = df)
    tibble(R2 = summary(m)$adj.r.squared, Pval = broom::glance(m)$p.value)
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
    }) %>%
    bind_rows()
  
  openxlsx::write.xlsx(stats_df,
                       paste0("Reco_Panel_Stats_R2_Pval_", model_label, ".xlsx"),
                       rowNames = FALSE)
  
  combined_preds <- left_join(combined_preds, stats_df, by = c("Plot_ID","Set","Model"))
  
  plot_scatter_by_set <- function(set_type) {
    df <- combined_preds %>% filter(Set == set_type, Model == model_label)
    annotation_df <- df %>%
      group_by(Plot_ID) %>% slice(1) %>%
      mutate(label = paste0("R² = ", round(R2, 2), "\nP = ", signif(Pval, 2)))
    ggplot(df, aes(x = Reco_obs, y = Reco_mod)) +
      geom_point(alpha = 0.6, color = "black", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, color = "gray40", linetype = "dashed") +
      geom_text(data = annotation_df,
                aes(x = Inf, y = Inf, label = label),
                inherit.aes = FALSE, size = 4.5, hjust = 1.1, vjust = 1.5) +
      facet_wrap(~ Plot_ID, ncol = 3, scales = "free") +
      labs(
        title = paste(set_type, "-", model_label, "- Observed vs Modeled Reco (by Plot)"),
        x = expression(Observed~Reco~(mu*mol~m^{-2}~s^{-1})),
        y = expression(Modeled~Reco~(mu*mol~m^{-2}~s^{-1}))
      ) +
      theme_bw(base_size = 12) +
      theme(strip.text = element_text(size = 11, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            axis.text  = element_text(size = 9),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank())
  }
  
  fig_cal  <- plot_scatter_by_set("Calibration")
  fig_eval <- plot_scatter_by_set("Evaluation")
  ggsave(paste0("Reco_Calibration_Scatterplots_byPlot_", model_label, ".png"),
         fig_cal,  width = 10, height = 7, dpi = 300)
  ggsave(paste0("Reco_Evaluation_Scatterplots_byPlot_",  model_label, ".png"),
         fig_eval, width = 10, height = 7, dpi = 300)
  
  # 6) Time series 3×2 (plot-level) — keep single PNG only (paginated PDF removed)
  data_ts <- Reco_predictions %>%
    mutate(Date = as.POSIXct(Date),
           Plot_ID = factor(Veg, levels = veg_levels_order)) %>%
    pivot_longer(cols = c(Reco_obs, Reco_mod),
                 names_to = "Type", values_to = "Reco") %>%
    filter(!is.na(Reco)) %>%
    mutate(Type = recode(Type,
                         Reco_obs = "Observed",
                         Reco_mod = "Modeled"))
  
  g_ts_all <- ggplot(data_ts, aes(x = Date, y = Reco, color = Type)) +
    geom_point(size = 1.4) +
    facet_wrap(~ Plot_ID, ncol = 3, nrow = 2, scales = "free_y") +
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
  
  # (Removed paginated PDF block)
  
  # 7) Build replicate-level predictions (B1..G3) from plot-level fit
  .parse_key_to_veg <- function(model_key) {
    mi <- strcapture("^([A-Z]+)_iter[0-9]+_(.+)$",
                     model_key,
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
        preds <- predict_manual(fit, df, temp_var = temp_var)
        out[[length(out) + 1]] <- data.frame(
          Veg      = veg,
          Rep      = rep_i,
          Plot_ID  = paste0(veg, rep_i),
          Model    = model_label,
          Date     = df$Date,
          Reco_obs = df$Reco,
          Reco_mod = preds,
          Residual = df$Reco - preds,
          Set      = set_name,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(out)) bind_rows(out) else
      data.frame(Veg=character(), Rep=character(), Plot_ID=character(),
                 Model=character(), Date=as.POSIXct(character()),
                 Reco_obs=numeric(), Reco_mod=numeric(), Residual=numeric(),
                 Set=character(), stringsAsFactors=FALSE)
  }
  
  cal_byrep  <- build_replicate_preds(Reco_cal,  "Calibration")
  eval_byrep <- build_replicate_preds(Reco_eval, "Evaluation")
  
  # 8) SCATTERPLOTS (BY REPLICATE)
  scatter_byrep <- function(df, title_txt) {
    stats_df <- df %>%
      group_by(Plot_ID) %>%
      summarise(
        R2   = { m <- try(lm(Reco_mod ~ Reco_obs, data = .), silent = TRUE);
        if (inherits(m, "try-error")) NA_real_ else summary(m)$adj.r.squared },
        Pval = { m <- try(lm(Reco_mod ~ Reco_obs, data = .), silent = TRUE);
        if (inherits(m, "try-error")) NA_real_ else broom::glance(m)$p.value },
        .groups = "drop"
      )
    df2 <- df %>%
      mutate(Plot_ID = factor(Plot_ID, levels = rep_levels_15))
    ann <- stats_df %>%
      mutate(label = paste0("R² = ", round(R2, 2), "\nP = ", signif(Pval, 2)))
    ggplot(df2, aes(x = Reco_obs, y = Reco_mod)) +
      geom_point(alpha = 0.6, size = 1.2, color = "black") +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray40") +
      geom_text(data = ann, aes(x = Inf, y = Inf, label = label),
                inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 4.2) +
      facet_wrap(~ Plot_ID, ncol = 3, scales = "free") +
      labs(title = title_txt,
           x = expression(Observed~Reco~(mu*mol~m^{-2}~s^{-1})),
           y = expression(Modeled~Reco~(mu*mol~m^{-2}~s^{-1}))) +
      theme_bw(base_size = 12) +
      theme(strip.text = element_text(size = 11, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank())
  }
  
  if (nrow(cal_byrep)) {
    p_cal_rep <- scatter_byrep(cal_byrep,
                               paste("Calibration —", model_label, "— Observed vs Modeled (by Replicate)"))
    ggsave(paste0("Reco_Calibration_Scatterplots_byReplicate_", model_label, ".png"),
           p_cal_rep, width = 10, height = 15, dpi = 300)
  }
  if (nrow(eval_byrep)) {
    p_eval_rep <- scatter_byrep(eval_byrep,
                                paste("Evaluation —", model_label, "— Observed vs Modeled (by Replicate)"))
    ggsave(paste0("Reco_Evaluation_Scatterplots_byReplicate_", model_label, ".png"),
           p_eval_rep, width = 10, height = 15, dpi = 300)
  }
  
  # 9) TIME SERIES 3×5 (BY REPLICATE)
  timeseries_byrep <- function(df, title_prefix) {
    df2 <- df %>%
      mutate(Plot_ID = factor(Plot_ID, levels = rep_levels_15)) %>%
      pivot_longer(c(Reco_obs, Reco_mod),
                   names_to = "Type", values_to = "Reco") %>%
      mutate(Type = recode(Type,
                           Reco_obs = "Observed",
                           Reco_mod = "Modeled"))
    ggplot(df2, aes(x = Date, y = Reco, color = Type)) +
      geom_point(size = 1.2, na.rm = TRUE) +
      facet_wrap(~ Plot_ID, ncol = 3, nrow = 5, scales = "free_y") +
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
    print(timeseries_byrep(cal_byrep,
                           "Observed vs Modeled – All 15 Replicates (Calibration)"))
    dev.off()
  }
  
  # 10) RMSE vs R² — Plot-level (REMOVED: no outputs saved)
  #    (Section intentionally deleted per request)
  
  # 11) Parameter histograms & summaries (a, b, c, d) — keep only per-plot PDF
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
                 names_to = "Parameter", values_to = "Estimate") %>%
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
    if (nrow(df) == 0) { warning(paste("No data for", param_name)); return(NULL) }
    ggplot(df, aes(x = SubsetType, y = Estimate)) +
      geom_col(color = "black", fill = NA, width = 0.7) +
      geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2) +
      geom_text(aes(label = Signif, y = Estimate + SE + 0.1), vjust = 0, size = 5, na.rm = TRUE) +
      facet_grid(cols = vars(Model), scales = "free_y") +
      labs(title = paste("Parameter", param_name, "by Plot (", model_label, ")", sep = ""),
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
    p_obj <- plot_param_hist(p); if (!is.null(p_obj)) print(p_obj)
  }
  dev.off()
  
  # (Averaged_Parameter_Histograms_* removed)
  # group summaries kept defined but no figure saved here
  long_params <- params %>%
    transmute(
      Group = Veg, Replicate = Veg, Model,
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
                 names_to = "Parameter", values_to = "Estimate") %>%
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
      .groups = "drop"
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
  
  plot_param_summary <- function(param_name) {
    df <- filter(summary_params, Parameter == param_name) %>% droplevels()
    ggplot(df, aes(x = Group, y = Mean_Estimate)) +
      geom_col(fill = "white", color = "black", width = 0.7) +
      geom_errorbar(aes(ymin = Mean_Estimate - SE_Estimate,
                        ymax = Mean_Estimate + SE_Estimate), width = 0.2) +
      geom_text(aes(label = Signif, y = Mean_Estimate + SE_Estimate + 0.05),
                size = 5, vjust = 0, na.rm = TRUE) +
      facet_grid(cols = vars(Model), scales = "free_y") +
      labs(title = paste("Parameter", param_name, "— group means (±SE),", model_label),
           x = "Vegetation Group", y = "Estimate ± SE") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold", size = 16),
            strip.text  = element_text(face = "bold", size = 12))
  }
  
  # (No ggsave/pdf for Averaged_Parameter_Histograms_*)
  
  # expose to global for final combined plot (combined panel removed, but assignments kept)
  assign(paste0("RESULTS_", model_label), results_all, envir = .GlobalEnv)
  assign(paste0("RECOCAL_",  model_label), Reco_cal,    envir = .GlobalEnv)
  assign(paste0("TEMPVAR_",  model_label), temp_var,     envir = .GlobalEnv)
  
  message("✅ Finished: ", model_label)
  invisible(list(
    results_all      = results_all,
    Reco_predictions = Reco_predictions,
    eval_results     = eval_results,
    cal_results      = cal_results,
    cal_byrep        = cal_byrep,
    eval_byrep       = eval_byrep,
    fitted_models    = fitted_models
  ))
}

# -----------------------------
# Run both models
# -----------------------------
out_Ts <- run_full_workflow("RE_5_Ts", "Model_Ts")
out_Ta <- run_full_workflow("RE_Ta",   "Model_Ta")

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

# ============================================================
# 2×2 PANEL — RMSE vs R² BY REPLICATE (Ta vs Ts)
# uses: out_Ta$cal_byrep, out_Ta$eval_byrep,
#       out_Ts$cal_byrep, out_Ts$eval_byrep
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(hydroGOF)
})

# ------------------------------------------------------------
# Helper: compute RMSE / adj.R² per replicate
# ------------------------------------------------------------
.compute_rep_metrics <- function(df, set_label = c("Calibration", "Evaluation")) {
  set_label <- match.arg(set_label)
  if (is.null(df) || !nrow(df)) return(data.frame())
  
  df %>%
    group_by(Veg, Rep, Plot_ID) %>%
    summarise(
      RMSE = hydroGOF::rmse(Reco_obs, Reco_mod),
      R2   = {
        ok <- is.finite(Reco_obs) & is.finite(Reco_mod)
        if (sum(ok) < 3) NA_real_ else summary(lm(Reco_obs[ok] ~ Reco_mod[ok]))$adj.r.squared
      },
      .groups = "drop"
    ) %>%
    mutate(Set = set_label)
}

# ------------------------------------------------------------
# Build metrics for each quadrant
# ------------------------------------------------------------
cal_rep_Ta  <- .compute_rep_metrics(out_Ta$cal_byrep,  "Calibration")
eval_rep_Ta <- .compute_rep_metrics(out_Ta$eval_byrep, "Evaluation")

cal_rep_Ts  <- .compute_rep_metrics(out_Ts$cal_byrep,  "Calibration")
eval_rep_Ts <- .compute_rep_metrics(out_Ts$eval_byrep, "Evaluation")

# ------------------------------------------------------------
# Helper: single panel RMSE vs R²
# ------------------------------------------------------------
.plot_rmse_r2_byrep <- function(df, title_txt) {
  if (!nrow(df)) {
    return(
      ggplot() + 
        ggtitle(paste(title_txt, "(no data)")) + 
        theme_void()
    )
  }
  
  # force order B, BG, BGS, GS, G
  df$Veg <- factor(df$Veg, levels = c("B","BG","BGS","GS","G"))
  df <- df %>% arrange(Veg, Rep)
  
  ggplot(df, aes(x = R2, y = RMSE, label = Plot_ID)) +
    geom_point(color = "black", size = 2.8) +
    geom_text_repel(size = 3, max.overlaps = 100, box.padding = 0.3, min.segment.length = 0) +
    scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0.02, 0.05))) +
    # unified y-scale so the 4 panels are comparable
    scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 0.5)) +
    labs(
      title = title_txt,
      x = expression(Adjusted~R^2),
      y = "RMSE (µmol m^{-2} s^{-1})"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------
# Build the 2×2 panel (Calibration/Evaluation × Ta/Ts)
# ------------------------------------------------------------
p_cal_ta  <- .plot_rmse_r2_byrep(cal_rep_Ta,  "Calibration — Air temperature (Ta)")
p_cal_ts  <- .plot_rmse_r2_byrep(cal_rep_Ts,  "Calibration — Soil temperature (Ts)")
p_eval_ta <- .plot_rmse_r2_byrep(eval_rep_Ta, "Evaluation — Air temperature (Ta)")
p_eval_ts <- .plot_rmse_r2_byrep(eval_rep_Ts, "Evaluation — Soil temperature (Ts)")

rmse_r2_byRep_Ta_vs_Ts_2x2 <- (p_cal_ta | p_cal_ts) /
  (p_eval_ta | p_eval_ts) +
  patchwork::plot_annotation(
    title = "RMSE vs adjusted R² by replicate for Ta and Ts models",
    subtitle = "Top: calibration • Bottom: evaluation • Left: Ta • Right: Ts",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
ggsave("RMSE_vs_R2_byRep_Ta_vs_Ts_2x2.png",
       rmse_r2_byRep_Ta_vs_Ts_2x2,
       width = 14, height = 10, dpi = 300)
ggsave("RMSE_vs_R2_byRep_Ta_vs_Ts_2x2.pdf",
       rmse_r2_byRep_Ta_vs_Ts_2x2,
       width = 14, height = 10)

cat("✅ Saved: RMSE_vs_R2_byRep_Ta_vs_Ts_2x2.(png|pdf)\n")
