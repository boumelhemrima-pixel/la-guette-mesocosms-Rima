# =============================
#      Load Required Libraries
# =============================
suppressPackageStartupMessages({
  library(readxl); library(openxlsx); library(dplyr)
  library(ggplot2); library(hydroGOF); library(ggrepel); library(tidyr)
})

# =============================
#           CONFIG
# =============================
setwd("C:/Users/rboumelh/Desktop/Mesocosms/Methane/Modeling CH4/3 temps_by treatment")

data_csv   <- "C:/Users/rboumelh/Desktop/Mesocosms/Modeling Reco/Full data_hourly averages.csv"
param_xlsx <- file.path(getwd(), "CH4_Model_Results_All_Model_Ts5.xlsx")  # ONE file with LAI params
Tref       <- 273.15 + 15  # Kelvin (15 °C)
setwd("C:/Users/rboumelh/Desktop/Mesocosms/Methane/Modeling CH4/3 temps_by treatment/gapfilling_Ts5")

# Canonical orders
group_levels <- c("B","BG","BGS","GS","G")
rep_levels   <- c("1","2","3")
plot_levels  <- as.vector(rbind(paste0("B", 1:3),
                                paste0("BG", 1:3),
                                paste0("BGS", 1:3),
                                paste0("GS", 1:3),
                                paste0("G", 1:3)))
groups <- group_levels
ids    <- plot_levels

# =============================
#         Helpers
# =============================
prefer_temp <- function(df, id) {
  Ts5 <- paste0("Ts5_", id)
  if (Ts5 %in% names(df)) Ts5 else "Tair_automatic"
}

# ER model (LAI only)
CH4_Model <- function(a, b, c, Ea, WTD, V2, Ts5, Tref) {
  ((a * (WTD / -0.35)) + (b * V2) + c) *
    exp((-Ea) * ((1 / (Ts5 + 273.15)) - (1 / Tref)))
}

# =============================
#  Load parameter file (LAI ONLY)
#  Expect columns: SubsetType, Model, a, b, c, Ea
#  Keep only rows with Model == "Model_Ts5"
# =============================
params_raw <- readxl::read_excel(param_xlsx)

needed <- c("Veg","Model","a","b","c","Ea")
miss   <- setdiff(needed, names(params_raw))
if (length(miss)) stop("Missing columns in parameter file: ", paste(miss, collapse=", "))

# Clean types
params_clean <- params_raw %>%
  dplyr::select(all_of(needed)) %>%
  dplyr::mutate(
    Veg   = as.character(Veg),
    Model = as.character(Model),
    a     = suppressWarnings(as.numeric(a)),
    b     = suppressWarnings(as.numeric(b)),
    c     = suppressWarnings(as.numeric(c)),
    Ea    = suppressWarnings(as.numeric(Ea))
  )

# Keep only the model we run
params_clean <- dplyr::filter(params_clean, Model == "Model_Ts5")
if (nrow(params_clean) == 0) stop("No rows with Model == 'Model_Ts5' found in parameter file.")

# ---- replicate Veg-level rows to B1..B3, BG1..BG3, ..., G3 ----
params_df <- params_clean %>%
  dplyr::filter(Veg %in% group_levels) %>%        # enforce canonical set/order
  tidyr::crossing(Rep = rep_levels) %>%           # 1,2,3 for each Veg
  dplyr::mutate(SubsetType = paste0(Veg, Rep)) %>%
  dplyr::select(SubsetType, Model, a, b, c, Ea)

# Build nested list used by run_model():
#   params_nested[["Model_Ts5"]][["B1"]] -> list(a=..., b=..., c=..., Ea=...)
params_by_model <- split(params_df, params_df$Model)
params_nested <- lapply(params_by_model, function(df) {
  setNames(
    lapply(seq_len(nrow(df)), function(i) as.list(df[i, c("a","b","c","Ea")])),
    df$SubsetType
  )
})

models_to_run <- names(params_nested)

# =============================
#    Run LAI model
# =============================
run_model <- function(model_name = "Model_Ts5") {
  message(sprintf(">>> Running %s", model_name))
  stopifnot(model_name %in% names(params_nested))
  
  # Load the environmental data once
  cbd_env <- read.csv(data_csv, sep = ",", stringsAsFactors = FALSE)
  cbd_env$Date_AVAL <- as.POSIXct(cbd_env$TIMESTAMP, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
  
  # LAI prefix only
  v2_prefix <- "LAI_"
  
  # Compute CH4_* columns per replicate id
  for (id in ids) {
    WTD_col <- paste0(id, "_WTD")
    V2_col  <- paste0(v2_prefix, id)         # "LAI_id"
    T_col   <- prefer_temp(cbd_env, id)
    CH4_col  <- paste0("CH4_", id)
    
    if (!all(c(WTD_col, V2_col, T_col) %in% names(cbd_env))) {
      warning(sprintf("[%s] missing one of [%s, %s, %s] → skipping %s",
                      model_name, WTD_col, V2_col, T_col, id))
      next
    }
    par_id <- params_nested[[model_name]][[id]]
    if (is.null(par_id)) {
      warning(sprintf("[%s] no parameters for %s in %s", model_name, id, basename(param_xlsx)))
      next
    }
    
    a  <- par_id$a; b <- par_id$b; c <- par_id$c; Ea <- par_id$Ea
    cbd_env[[CH4_col]] <- CH4_Model(
      a, b, c, Ea,
      WTD = cbd_env[[WTD_col]],
      V2  = cbd_env[[V2_col]],
      Ts5 = cbd_env[[T_col]],
      Tref= Tref
    )
  }
  
  # Keep negatives (NO clamping)
  CH4_cols <- grep("^CH4_", names(cbd_env), value = TRUE)
  
  # -----------------------------
  # Save hourly with CH4 columns
  # -----------------------------
  out_tag <- model_name
  dir.create(file.path("gapfilling", out_tag), recursive = TRUE, showWarnings = FALSE)
  write.csv(cbd_env, sprintf("gapfilling/%s/cbd_CH4_gap_filled_%s.csv", out_tag, out_tag), row.names = FALSE)
  
  # -----------------------------
  # Daily aggregation
  # -----------------------------
  cbd_env$Date_DAY <- as.Date(cbd_env$Date_AVAL, tz = "Europe/Paris")
  Nm <- aggregate(cbd_env[CH4_cols], by = list(Date_AVAL = cbd_env$Date_DAY), FUN = sum, na.rm = TRUE)
  
  # Daily means & std by vegetation group
  combo <- list(
    B   = c("CH4_B1","CH4_B2","CH4_B3"),
    BG  = c("CH4_BG1","CH4_BG2","CH4_BG3"),
    BGS = c("CH4_BGS1","CH4_BGS2","CH4_BGS3"),
    GS  = c("CH4_GS1","CH4_GS2","CH4_GS3"),
    G   = c("CH4_G1","CH4_G2","CH4_G3")
  )
  for (g in names(combo)) {
    Nm[[paste0("Mean_CH4_", g)]] <- rowMeans(Nm[, combo[[g]], drop=FALSE], na.rm = TRUE)
    Nm[[paste0("Std_CH4_",  g)]] <- apply(Nm[, combo[[g]], drop=FALSE], 1, sd, na.rm = TRUE)
  }
  
  # -----------------------------
  # Time windows
  # -----------------------------
  start_2023 <- as.Date("2023-02-01")
  end_2024   <- as.Date("2024-02-01")
  end_2025   <- as.Date("2025-02-01")
  
  Nmm <- subset(Nm, Date_AVAL >= start_2023 & Date_AVAL < end_2025)
  
  y_lim <- c(-0.01, 0.2)
  xlim_dates <- as.Date(c("2023-02-01", "2025-02-01"))
  x_ticks    <- as.Date(c("2023-02-01", "2024-02-01", "2025-02-01"))
  x_labels   <- c("2023","2024","2025")
  
  # =========================================================
  # FIGURE: Hourly replicate fluxes (3 × 5 grid) — MODELED only
  # =========================================================
  rep_hourly <- do.call(rbind, lapply(groups, function(g) {
    do.call(rbind, lapply(1:3, function(rp) {
      cn <- paste0("CH4_", g, rp)
      if (!cn %in% names(cbd_env)) return(NULL)
      data.frame(
        Date    = cbd_env$Date_AVAL,
        Group   = g,
        Rep     = factor(rp, levels = rep_levels),
        Plot_ID = paste0(g, rp),
        CH4    = cbd_env[[cn]],  # hourly µmol m^-2 s^-1
        stringsAsFactors = FALSE
      )
    }))
  }))
  rep_hourly <- rep_hourly[!is.na(rep_hourly$CH4) & !is.na(rep_hourly$Date), , drop = FALSE]
  rep_hourly$Group   <- factor(rep_hourly$Group,   levels = group_levels)
  rep_hourly$Plot_ID <- factor(rep_hourly$Plot_ID, levels = plot_levels)
  rep_hourly$Rep     <- factor(rep_hourly$Rep,     levels = rep_levels)
  
  if (nrow(rep_hourly) > 0) {
    p_rep_hourly <- ggplot(rep_hourly, aes(Date, CH4)) +
      geom_point(size = 0.3, alpha = 0.6) +
      facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
      coord_cartesian(ylim = y_lim) +
      labs(title = "Hourly CH4 per replicate – Model_Ts5",
           x = "Date", y = expression(R[eco]~(mu*mol~m^{-2}~s^-1))) +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(face = "bold", size = 8))
    
    ggsave(sprintf("gapfilling/%s/Replicate_CH4_all_years_HOURLY_%s.png", out_tag, out_tag),
           plot = p_rep_hourly, width = 12, height = 15, dpi = 300)
  }
  
  # =========================================================
  # FIGURE: Hourly CH4 split into DAY vs NIGHT (PAR=0 → Night)
  # =========================================================
  if ("PAR" %in% names(cbd_env)) {
    rep_hourly_dn <- rep_hourly %>%
      left_join(cbd_env %>% select(Date_AVAL, PAR), by = c("Date" = "Date_AVAL")) %>%
      mutate(DayNight = ifelse(PAR == 0, "Night", "Day"))
    
    if (nrow(rep_hourly_dn) > 0) {
      p_rep_dn <- ggplot(rep_hourly_dn, aes(Date, CH4, color = DayNight)) +
        geom_point(size = 0.3, alpha = 0.6) +
        facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
        coord_cartesian(ylim = y_lim) +
        scale_color_manual(values = c("Day" = "orange", "Night" = "blue")) +
        labs(title = "Hourly CH4 per replicate – Model_Ts5 (Day vs Night by PAR)",
             x = "Date", y = expression(R[eco]~(mu*mol~m^{-2}~s^-1)), color = "") +
        theme_bw(base_size = 10) +
        theme(strip.text = element_text(face = "bold", size = 8),
              legend.position = "top")
      
      ggsave(sprintf("gapfilling/%s/Replicate_CH4_all_years_HOURLY_DayNight_%s.png", 
                     out_tag, out_tag),
             plot = p_rep_dn, width = 12, height = 15, dpi = 300)
    }
  } else {
    warning("PAR column not found in cbd_env → cannot split Day/Night")
  }
  
  # =========================================================
  # FIGURE: Hourly replicate CH4 – DAY only
  # =========================================================
  if ("PAR" %in% names(cbd_env) && nrow(rep_hourly) > 0) {
    rep_day <- rep_hourly %>%
      filter(cbd_env$PAR[match(Date, cbd_env$Date_AVAL)] > 0)
    
    if (nrow(rep_day) > 0) {
      p_rep_day <- ggplot(rep_day, aes(Date, CH4)) +
        geom_point(size = 0.3, alpha = 0.6, color = "orange") +
        facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
        coord_cartesian(ylim = y_lim) +
        labs(title = "Hourly CH4 per replicate – DAY only (Model_Ts5)",
             x = "Date", y = expression(R[eco]~(mu*mol~m^{-2}~s^-1))) +
        theme_bw(base_size = 10) +
        theme(strip.text = element_text(face = "bold", size = 8))
      
      ggsave(sprintf("gapfilling/%s/Replicate_CH4_all_years_HOURLY_DAY_%s.png", out_tag, out_tag),
             plot = p_rep_day, width = 12, height = 15, dpi = 300)
    }
  }
  
  # =========================================================
  # FIGURE: Hourly replicate CH4 – NIGHT only
  # =========================================================
  if ("PAR" %in% names(cbd_env) && nrow(rep_hourly) > 0) {
    rep_night <- rep_hourly %>%
      filter(cbd_env$PAR[match(Date, cbd_env$Date_AVAL)] == 0)
    
    if (nrow(rep_night) > 0) {
      p_rep_night <- ggplot(rep_night, aes(Date, CH4)) +
        geom_point(size = 0.3, alpha = 0.6, color = "blue") +
        facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
        coord_cartesian(ylim = y_lim) +
        labs(title = "Hourly CH4 per replicate – NIGHT only (Model_Ts5)",
             x = "Date", y = expression(R[eco]~(mu*mol~m^{-2}~s^-1))) +
        theme_bw(base_size = 10) +
        theme(strip.text = element_text(face = "bold", size = 8))
      
      ggsave(sprintf("gapfilling/%s/Replicate_CH4_all_years_HOURLY_NIGHT_%s.png", out_tag, out_tag),
             plot = p_rep_night, width = 12, height = 15, dpi = 300)
    }
  }
  
  # =========================================================
  # FIGURE: Hourly observed vs modeled (3 × 5)
  # =========================================================
  obs_csv <- "C:/Users/rboumelh/Desktop/Mesocosms/Fluxes/mesocosm_data1.csv"
  obs_df  <- tryCatch(
    read.csv(obs_csv, sep = ",", stringsAsFactors = FALSE),
    error = function(e) { warning("Could not read obs file: ", conditionMessage(e)); NULL }
  )
  obs_long <- NULL
  if (!is.null(obs_df)) {
    # Date parsing
    if ("Date" %in% names(obs_df)) {
      obs_df$Date <- as.POSIXct(obs_df$Date, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
    } else if ("TIMESTAMP" %in% names(obs_df)) {
      obs_df$Date <- as.POSIXct(obs_df$TIMESTAMP, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
    } else if ("TIMESTAMP.Reco" %in% names(obs_df)) {
      obs_df$Date <- as.POSIXct(obs_df$TIMESTAMP.Reco, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
    } else {
      warning("No recognizable Date column in obs_df; skipping overlay.")
      obs_df <- NULL
    }
  }
  if (!is.null(obs_df)) {
    # Wide schema preferred: CH4_B1, CH4_BG2, ...
    wide_cols <- paste0("CH4_", plot_levels)
    if (any(wide_cols %in% names(obs_df))) {
      obs_long <- do.call(rbind, lapply(groups, function(g) {
        do.call(rbind, lapply(1:3, function(rp) {
          cn <- paste0("CH4_", g, rp)
          if (!cn %in% names(obs_df)) return(NULL)
          data.frame(
            Date     = obs_df$Date,
            Group    = g,
            Rep      = factor(rp, levels = rep_levels),
            Plot_ID  = paste0(g, rp),
            CH4_obs = obs_df[[cn]],
            stringsAsFactors = FALSE
          )
        }))
      }))
    } else {
      # Long schema: Plot/Mesocosm + Rep/Replicates + CH4
      plot_col <- if ("Plot" %in% names(obs_df)) "Plot" else if ("Mesocosm" %in% names(obs_df)) "Mesocosm" else NA
      rep_col  <- if ("Rep"  %in% names(obs_df)) "Rep"  else if ("Replicates" %in% names(obs_df)) "Replicates" else NA
      if (!is.na(plot_col) && !is.na(rep_col) && "CH4" %in% names(obs_df)) {
        tmp <- obs_df[, c("Date", plot_col, rep_col, "CH4")]
        names(tmp) <- c("Date", "Group", "Rep", "CH4_obs")
        tmp$Group   <- as.character(tmp$Group)
        tmp$Rep     <- as.character(tmp$Rep)
        tmp$Plot_ID <- paste0(tmp$Group, tmp$Rep)
        tmp <- tmp[tmp$Plot_ID %in% plot_levels, , drop = FALSE]
        obs_long <- tmp
      }
    }
  }
  if (!is.null(obs_long) && nrow(obs_long) > 0 && nrow(rep_hourly) > 0) {
    obs_long$Group   <- factor(obs_long$Group,   levels = group_levels)
    obs_long$Rep     <- factor(obs_long$Rep,     levels = rep_levels)
    obs_long$Plot_ID <- factor(obs_long$Plot_ID, levels = plot_levels)
    
    rep_hourly_merged <- merge(rep_hourly, obs_long,
                               by = c("Date","Group","Rep","Plot_ID"),
                               all.x = TRUE)
    
    p_rep_hourly_ovl <- ggplot(rep_hourly_merged, aes(Date)) +
      # draw modeled first (goes to the "background")
      geom_point(aes(y = CH4, color = "Modeled"),
                 size = 0.3, alpha = 0.5, na.rm = TRUE) +
      # draw observed second (always on top)
      geom_point(aes(y = CH4_obs, color = "Observed"),
                 size = 1, alpha = 1, na.rm = TRUE) +
      facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
      coord_cartesian(ylim = y_lim) +
      scale_color_manual(values = c("Observed" = "red", "Modeled" = "black")) +
      labs(title = "Hourly observed vs modeled CH4 per replicate – Model_Ts5",
           x = "Date", y = expression(R[eco]~(mu*mol~m^{-2}~s^-1)), color = "") +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(face = "bold", size = 8),
            legend.position = "top")
    
    ggsave(sprintf("gapfilling/%s/Replicate_CH4_all_years_HOURLY_vsOBS_%s.png", out_tag, out_tag),
           plot = p_rep_hourly_ovl, width = 12, height = 15, dpi = 300)
  } else {
    warning("Observed overlay skipped: no observed data found or schema not recognized.")
  }
  
  # =========================================================
  # FIGURE: ONLY corresponding timestamps (Observed & Modeled)
  # =========================================================
  if (exists("rep_hourly_merged")) {
    corr_df <- rep_hourly_merged %>%
      dplyr::filter(is.finite(CH4), is.finite(CH4_obs))
    
    if (nrow(corr_df) > 0) {
      p_corr <- ggplot(corr_df, aes(Date)) +
        geom_point(aes(y = CH4_obs, color = "Observed"),
                   size = 2, alpha = 1, na.rm = TRUE) +
        geom_point(aes(y = CH4,     color = "Modeled"),
                   size = 2, alpha = 1, na.rm = TRUE) +
        facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
        coord_cartesian(ylim = y_lim) +
        scale_color_manual(values = c("Observed" = "red", "Modeled" = "black")) +
        labs(title = "Hourly CH4 — corresponding observed & gap-filled only (Model_Ts5)",
             x = "Date", y = expression(R[eco]~(mu*mol~m^{-2}~s^{-1})), color = "") +
        theme_bw(base_size = 10) +
        theme(strip.text = element_text(face = "bold", size = 8),
              legend.position = "top")
      
      ggsave(sprintf("gapfilling/%s/Replicate_CH4_HOURLY_vsOBS_CORRESP_%s.png", out_tag, out_tag),
             plot = p_corr, width = 12, height = 15, dpi = 300)
    } else {
      warning("No corresponding timestamps found where both observed and modeled are present.")
    }
  }
  # =========================================================
  
  # -----------------------------
  # Yearly sums per replicate + Excel
  # -----------------------------
  summarize_year <- function(df_year, tag) {
    num_cols <- vapply(df_year, is.numeric, TRUE)
    df_year[num_cols] <- lapply(df_year[num_cols], function(x) { x[!is.finite(x)] <- NA; x })
    write.csv(df_year, sprintf("gapfilling/%s/cbd_env%s_%s.csv", out_tag, tag, out_tag), row.names = FALSE)
    
    sums <- lapply(unlist(combo, use.names = FALSE), function(col) {
      sum(df_year[[col]], na.rm=TRUE) * (3600 * 12 / 1e6)  # keep your unit conversion
    })
    names(sums) <- unlist(combo, use.names = TRUE)
    
    m <- data.frame(
      Category = names(combo),
      Mean = sapply(names(combo), function(g) mean(unlist(sums[combo[[g]]]))),
      SD   = sapply(names(combo), function(g)   sd(unlist(sums[combo[[g]]])))
    )
    list(meansd = m, sums = sums)
  }
  
  year_2023 <- subset(Nm, Date_AVAL >= start_2023 & Date_AVAL < end_2024)
  year_2024 <- subset(Nm, Date_AVAL >= end_2024   & Date_AVAL < end_2025)
  
  res23 <- summarize_year(year_2023, "2023")
  res24 <- summarize_year(year_2024, "2024")
  
  wb <- createWorkbook()
  addWorksheet(wb, "2023"); writeData(wb, "2023", res23$meansd)
  addWorksheet(wb, "2024"); writeData(wb, "2024", res24$meansd)
  addWorksheet(wb, "All_Replicates")
  
  all_repl_df <- function(sums_list, yr) {
    data.frame(
      Year      = yr,
      Category  = rep(names(combo), each = 3),
      Replicate = rep(1:3, times = length(combo)),
      Value     = unlist(lapply(combo, function(cols) unname(unlist(sums_list[cols]))))
    )
  }
  all_rep <- rbind(all_repl_df(res23$sums, "2023"),
                   all_repl_df(res24$sums, "2024"))
  writeData(wb, "All_Replicates", all_rep)
  saveWorkbook(wb, sprintf("gapfilling/%s/CH4_Results_%s.xlsx", out_tag, out_tag), overwrite = TRUE)
  
  invisible(list(Nm = Nm, Nmm = Nmm, rep_hourly = rep_hourly))
}

# =============================
#         RUN (LAI ONLY)
# =============================
dir.create(file.path(getwd(), "gapfilling", "Model_Ts5"), recursive = TRUE, showWarnings = FALSE)
results <- list(Model_Ts5 = tryCatch(run_model("Model_Ts5"),
                                     error = function(e) { message("!! Model_Ts5 failed: ", conditionMessage(e)); NULL }))

