# =============================
#      Load Required Libraries
# =============================
suppressPackageStartupMessages({
  library(readxl); library(openxlsx); library(dplyr)
  library(ggplot2); library(tidyr)
})

# =============================
#           CONFIG
# =============================
setwd("C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by 2nd method/by treatment")

data_csv   <- "C:/Users/rboumelh/Desktop/Mesocosms/Modeling Reco/Full data_hourly averages.csv"
param_xlsx <- file.path(getwd(), "GPP_Model_Results_All_LAI_Groups_withReplicates.xlsx")  # a,b,Topt per B1..G3

# Canonical orders
group_levels <- c("B","BG","BGS","GS","G")
rep_levels   <- c("1","2","3")
plot_levels  <- as.vector(rbind(paste0("B", 1:3),
                                paste0("BG", 1:3),
                                paste0("BGS", 1:3),
                                paste0("GS", 1:3),
                                paste0("G", 1:3)))

# Vegetation groups / replicate IDs (forced to canonical order)
groups <- group_levels
ids    <- plot_levels

# =============================
#         Helpers
# =============================
prefer_temp <- function(df, id) {
  ta <- paste0("Ta_", id)
  if (ta %in% names(df)) ta else "Tair_automatic"
}

# Bounded temperature response (Tmin=0, Tmax=40)
fT_bounded <- function(Ta, Topt, Tmin = 0, Tmax = 40) {
  num   <- (Ta - Tmin) * (Ta - Tmax)
  denom <- (num - (Ta - Topt)^2)
  out <- num / denom
  out[!is.finite(out)] <- NA_real_
  out
}

# GPP model for either VI or LAI (V2 is VI_* or LAI_* column)
gpp_Model <- function(a, b, Topt, K, V2, Ta, PAR) {
  fT <- fT_bounded(Ta, Topt)
  (((a * V2 + b) * fT) * PAR) / (K + PAR)
}

# =============================
#  Load SINGLE parameter file
#  Expect columns: SubsetType (e.g., B1), Model (Model_VI/Model_LAI), a, b, Topt
# =============================
# =============================
#  Load SINGLE parameter file (Veg-level) and replicate to B1..G3
#  Expected columns in Excel: Veg, Model, a, b, Topt, K (others ignored)
# =============================

params_raw <- readxl::read_excel(param_xlsx)

needed0 <- c("Veg","Model","a","b","Topt","K")
miss0   <- setdiff(needed0, names(params_raw))
if (length(miss0)) stop("Missing columns in parameter file: ", paste(miss0, collapse=", "))

# Keep only the columns we need, ensure proper types
params_clean <- params_raw %>%
  dplyr::select(all_of(needed0)) %>%
  dplyr::mutate(
    Veg   = as.character(Veg),
    Model = as.character(Model),
    a     = suppressWarnings(as.numeric(a)),
    b     = suppressWarnings(as.numeric(b)),
    Topt  = suppressWarnings(as.numeric(Topt)),
    K     = suppressWarnings(as.numeric(K))
  )

# Normalize model names so your code keeps working
# (screenshot shows "LAI"; we convert to "Model_LAI")
params_clean$Model <- dplyr::recode(
  params_clean$Model,
  "LAI" = "Model_LAI",
  "VI"  = "Model_VI",
  .default = params_clean$Model
)

# Keep only models we actually run (LAI here)
params_clean <- dplyr::filter(params_clean, Model %in% c("Model_LAI","Model_VI"))

# --- NEW: replicate Veg-level rows to B1..B3 etc. ---
params_df <- params_clean %>%
  dplyr::filter(Veg %in% group_levels) %>%                # enforce canonical order
  tidyr::crossing(Rep = rep_levels) %>%                   # make 1,2,3 for each Veg
  dplyr::mutate(SubsetType = paste0(Veg, Rep)) %>%        # B1, B2, ...
  dplyr::select(SubsetType, Model, a, b, Topt, K)

# Build nested parameter list used by the model runner
params_by_model <- split(params_df, params_df$Model)
params_nested <- lapply(params_by_model, function(df) {
  setNames(
    lapply(seq_len(nrow(df)), function(i) as.list(df[i, c("a","b","Topt","K")])),
    df$SubsetType
  )
})

models_to_run <- names(params_nested)
if (!length(models_to_run)) stop("No usable rows found in parameter file.")

# =============================
#    Run ONE GPP model
# =============================
run_model <- function(model_name) {
  message(sprintf(">>> Running %s", model_name))
  stopifnot(model_name %in% names(params_nested))
  
  # Load environmental data
  cbd_env <- read.csv(data_csv, sep = ",", stringsAsFactors = FALSE)
  cbd_env$Date_AVAL <- as.POSIXct(cbd_env$TIMESTAMP, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
  
  # Choose VI/LAI prefix by model
  v2_prefix <- if (grepl("LAI$", model_name)) "LAI_" else "VI_"
  
  # Ensure PAR is numeric
  if (!"PAR" %in% names(cbd_env)) stop("PAR column not found in environmental data.")
  cbd_env$PAR <- suppressWarnings(as.numeric(cbd_env$PAR))
  
  # Compute GPP_* columns per replicate id
  for (id in ids) {
    V2_col <- paste0(v2_prefix, id)      # "VI_B1" or "LAI_B1", ...
    T_col  <- prefer_temp(cbd_env, id)   # "Ta_B1" or "Tair_automatic"
    GPP_col <- paste0("GPP_", id)
    
    if (!all(c(V2_col, T_col) %in% names(cbd_env))) {
      warning(sprintf("[%s] missing one of [%s, %s] → skipping %s",
                      model_name, V2_col, T_col, id))
      next
    }
    par_id <- params_nested[[model_name]][[id]]
    if (is.null(par_id)) {
      warning(sprintf("[%s] no parameters for %s in %s", model_name, id, basename(param_xlsx)))
      next
    }
    
    a  <- par_id$a; b <- par_id$b; Topt <- par_id$Topt; K <- par_id$K
    cbd_env[[GPP_col]] <- gpp_Model(
      a, b, Topt, K,
      V2  = suppressWarnings(as.numeric(cbd_env[[V2_col]])),
      Ta  = suppressWarnings(as.numeric(cbd_env[[T_col]])),
      PAR = cbd_env$PAR
    )
  }
  
  # Columns created
  gpp_cols <- grep("^GPP_", names(cbd_env), value = TRUE)
  
  # -----------------------------
  # Save hourly with GPP columns
  # -----------------------------
  out_tag <- model_name
  dir.create(file.path("gapfilling", out_tag), recursive = TRUE, showWarnings = FALSE)
  write.csv(cbd_env, sprintf("gapfilling/%s/cbd_GPP_gap_filled_%s.csv", out_tag, out_tag), row.names = FALSE)
  
  # -----------------------------
  # Daily aggregation (sum) then convert to mean hourly (/24)
  # -----------------------------
  cbd_env$Date_DAY <- as.Date(cbd_env$Date_AVAL, tz = "Europe/Paris")
  Nm <- aggregate(cbd_env[gpp_cols], by = list(Date_AVAL = cbd_env$Date_DAY), FUN = sum, na.rm = TRUE)
  
  # Daily means & std by vegetation group (across 3 reps)
  combo <- list(
    B   = c("GPP_B1","GPP_B2","GPP_B3"),
    BG  = c("GPP_BG1","GPP_BG2","GPP_BG3"),
    BGS = c("GPP_BGS1","GPP_BGS2","GPP_BGS3"),
    GS  = c("GPP_GS1","GPP_GS2","GPP_GS3"),
    G   = c("GPP_G1","GPP_G2","GPP_G3")
  )
  for (g in names(combo)) {
    Nm[[paste0("Mean_GPP_", g)]] <- rowMeans(Nm[, combo[[g]], drop=FALSE], na.rm = TRUE)
    Nm[[paste0("Std_GPP_",  g)]] <- apply(Nm[, combo[[g]], drop=FALSE], 1, sd, na.rm = TRUE)
  }
  
  # -----------------------------
  # Time windows
  # -----------------------------
  start_2023 <- as.Date("2023-02-01")
  end_2024   <- as.Date("2024-02-01")
  end_2025   <- as.Date("2025-02-01")
  
  Nmm <- subset(Nm, Date_AVAL >= start_2023 & Date_AVAL < end_2025)

  # =========================================================
  # FIGURE: Hourly GPP per replicate — MODELED only
  # =========================================================
  rep_hourly <- do.call(rbind, lapply(groups, function(g) {
    do.call(rbind, lapply(1:3, function(rp) {
      cn <- paste0("GPP_", g, rp)
      if (!cn %in% names(cbd_env)) return(NULL)
      data.frame(
        Date    = cbd_env$Date_AVAL,
        Group   = g,
        Rep     = factor(rp, levels = rep_levels),
        Plot_ID = paste0(g, rp),
        GPP     = cbd_env[[cn]],  # hourly modeled GPP (µmol CO2 m^-2 s^-1)
        stringsAsFactors = FALSE
      )
    }))
  }))
  rep_hourly <- rep_hourly[!is.na(rep_hourly$GPP) & !is.na(rep_hourly$Date), , drop = FALSE]
  rep_hourly$Group   <- factor(rep_hourly$Group,   levels = group_levels)
  rep_hourly$Plot_ID <- factor(rep_hourly$Plot_ID, levels = plot_levels)
  rep_hourly$Rep     <- factor(rep_hourly$Rep,     levels = rep_levels)
  
  y_lim_hourly <- c(-35, 0)  # same sign convention as your daily mean plot
  
  if (nrow(rep_hourly) > 0) {
    p_rep_hourly <- ggplot(rep_hourly, aes(Date, GPP)) +
      geom_point(size = 0.3, alpha = 0.6) +
      facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
      coord_cartesian(ylim = y_lim_hourly) +
      labs(title = paste("Hourly GPP per replicate –", model_name),
           x = "Date", y = expression("GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")")) +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(face = "bold", size = 8))
    
    ggsave(sprintf("gapfilling/%s/Replicate_GPP_all_years_HOURLY_%s.png", out_tag, out_tag),
           plot = p_rep_hourly, width = 12, height = 15, dpi = 300)
  }
  
  # =========================================================
  # FIGURE: Hourly observed vs modeled (3 × 5)
  #   - Observed GPP auto-detected (wide or long schema)
  # =========================================================
  obs_csv_candidates <- c(
    "C:/Users/rboumelh/Desktop/Mesocosms/Fluxes/mesocosm_data1.csv",
    "C:/Users/rboumelh/Desktop/Mesocosms/Modeling Reco/mesocosm_data1.csv"
  )
  obs_csv <- obs_csv_candidates[which(file.exists(obs_csv_candidates))][1]
  if (is.na(obs_csv)) {
    warning("No observed GPP file found in candidates; skipping overlay.")
  }
  
  obs_df  <- if (!is.na(obs_csv)) {
    tryCatch(read.csv(obs_csv, sep = ",", stringsAsFactors = FALSE),
             error = function(e) { warning("Could not read obs file: ", conditionMessage(e)); NULL })
  } else NULL
  
  # Try to create a proper Date column
  if (!is.null(obs_df)) {
    if ("Date" %in% names(obs_df)) {
      obs_df$Date <- as.POSIXct(obs_df$Date, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
    } else if ("TIMESTAMP" %in% names(obs_df)) {
      obs_df$Date <- as.POSIXct(obs_df$TIMESTAMP, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
    } else if ("TIMESTAMP.Reco" %in% names(obs_df)) {
      obs_df$Date <- as.POSIXct(obs_df$TIMESTAMP.Reco, format = "%m/%d/%Y %H:%M", tz = "Europe/Paris")
    } else {
      warning("No recognizable Date column in observed file; skipping overlay.")
      obs_df <- NULL
    }
  }
  
  obs_long <- NULL
  if (!is.null(obs_df)) {
    # (A) Wide schema: GPP_B1, GPP_BG2, ...
    wide_cols <- paste0("GPP_", plot_levels)
    if (any(wide_cols %in% names(obs_df))) {
      obs_long <- do.call(rbind, lapply(groups, function(g) {
        do.call(rbind, lapply(1:3, function(rp) {
          cn <- paste0("GPP_", g, rp)
          if (!cn %in% names(obs_df)) return(NULL)
          data.frame(
            Date     = obs_df$Date,
            Group    = g,
            Rep      = factor(rp, levels = rep_levels),
            Plot_ID  = paste0(g, rp),
            GPP_obs  = obs_df[[cn]],
            stringsAsFactors = FALSE
          )
        }))
      }))
    } else {
      # (B) Long schema: Plot/Mesocosm + Rep/Replicates + GPP
      plot_col <- if ("Plot" %in% names(obs_df)) "Plot" else if ("Mesocosm" %in% names(obs_df)) "Mesocosm" else NA
      rep_col  <- if ("Rep"  %in% names(obs_df)) "Rep"  else if ("Replicates" %in% names(obs_df)) "Replicates" else NA
      val_col  <- if ("GPP"  %in% names(obs_df)) "GPP"  else NA
      if (!is.na(plot_col) && !is.na(rep_col) && !is.na(val_col)) {
        tmp <- obs_df[, c("Date", plot_col, rep_col, val_col)]
        names(tmp) <- c("Date", "Group", "Rep", "GPP_obs")
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
    
    p_hourly_ovl <- ggplot(rep_hourly_merged, aes(Date)) +
      # Modeled first (background)
      geom_point(aes(y = GPP,     color = "Modeled"),
                 size = 0.3, alpha = 0.55, na.rm = TRUE) +
      # Observed on top
      geom_point(aes(y = GPP_obs, color = "Observed"),
                 size = 1, alpha = 1, na.rm = TRUE) +
      facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
      coord_cartesian(ylim = y_lim_hourly) +
      scale_color_manual(values = c("Observed" = "red", "Modeled" = "black")) +
      labs(title = paste("Hourly observed vs modeled GPP per replicate –", model_name),
           x = "Date", y = expression("GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"), color = "") +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(face = "bold", size = 8),
            legend.position = "top")
    
    ggsave(sprintf("gapfilling/%s/Replicate_GPP_all_years_HOURLY_vsOBS_%s.png", out_tag, out_tag),
           plot = p_hourly_ovl, width = 12, height = 15, dpi = 300)
    
    # =========================================================
    # FIGURE: ONLY corresponding timestamps (Observed & Modeled)
    # =========================================================
    corr_df <- rep_hourly_merged %>%
      dplyr::filter(is.finite(GPP), is.finite(GPP_obs))
    
    if (nrow(corr_df) > 0) {
      p_corr <- ggplot(corr_df, aes(Date)) +
        geom_point(aes(y = GPP_obs, color = "Observed"),
                   size = 2, alpha = 0.75, na.rm = TRUE) +
        geom_point(aes(y = GPP,     color = "Modeled"),
                   size = 2, alpha = 0.6, na.rm = TRUE) +
        facet_grid(rows = vars(Group), cols = vars(Rep), scales = "free_y") +
        coord_cartesian(ylim = y_lim_hourly) +
        scale_color_manual(values = c("Observed" = "red", "Modeled" = "black")) +
        labs(title = paste("Hourly GPP — corresponding observed & gap-filled only (", model_name, ")", sep = ""),
             x = "Date", y = expression("GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"), color = "") +
        theme_bw(base_size = 10) +
        theme(strip.text = element_text(face = "bold", size = 8),
              legend.position = "top")
      
      ggsave(sprintf("gapfilling/%s/Replicate_GPP_HOURLY_vsOBS_CORRESP_%s.png", out_tag, out_tag),
             plot = p_corr, width = 12, height = 15, dpi = 300)
    } else {
      warning("No corresponding timestamps where both observed and modeled GPP are present.")
    }
  } else {
    warning("Observed overlay skipped: no observed GPP data found or schema not recognized.")
  }
  # =========================================================
  
  # -----------------------------
  # Yearly sums per replicate + Excel
  # -----------------------------
  summarize_year <- function(df_year, tag) {
    num_cols <- vapply(df_year, is.numeric, TRUE)
    df_year[num_cols] <- lapply(df_year[num_cols], function(x) { x[!is.finite(x)] <- NA; x })
    write.csv(df_year, sprintf("gapfilling/%s/cbd_env_%s_%s.csv", out_tag, tag, out_tag), row.names = FALSE)
    
    sums <- lapply(unlist(combo, use.names = FALSE), function(col) {
      sum(df_year[[col]], na.rm=TRUE) * (3600 * 12 / 1e6)  # (keep your unit conversion)
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
  saveWorkbook(wb, sprintf("gapfilling/%s/GPP_Results_%s.xlsx", out_tag, out_tag), overwrite = TRUE)
  
  invisible(list(Nm = Nm, Nmm = Nmm))
}

# =============================
#         RUN MODELS
# =============================
for (mn in models_to_run) {
  dir.create(file.path(getwd(), "gapfilling", mn), recursive = TRUE, showWarnings = FALSE)
}
results <- setNames(vector("list", length(models_to_run)), models_to_run)
for (mn in models_to_run) {
  results[[mn]] <- tryCatch(run_model(mn),
                            error = function(e) { message("!! ", mn, " failed: ", conditionMessage(e)); NULL })
}

