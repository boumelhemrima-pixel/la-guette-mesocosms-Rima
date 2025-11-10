# =============================
#      Load Required Libraries
# =============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(hydroGOF)
  library(scales)
  library(cowplot)   # for get_legend / plot_grid
  library(grid)
  library(openxlsx)  # <-- added
})

# =============================
#           CONFIG
# =============================
setwd("C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by campaign/by replicate/GPP_by treatment/for git")

in_csv  <- "mesocosm_data1_with_GPPmax_predictions.csv"
out_rep <- "Obs_vs_Modeled_GPP_byReplicate_VI_LAI_3x5.png"
out_veg <- "Obs_vs_Modeled_GPP_byVeg_VI_LAI_3x2.png"

# master toggle (keep code, but disable VI if FALSE)
USE_VI <- FALSE

# Fixed facet orders
meso_order <- c("B1","B2","B3",
                "BG1","BG2","BG3",
                "BGS1","BGS2","BGS3",
                "GS1","GS2","GS3",
                "G1","G2","G3")
veg_levels <- c("B","BG","BGS","GS","G")

# Per-replicate K values (your list)
k_values <- c(
  "B1" = 407.0161559, "B2" = 396.1684879, "B3" = 744.6314954,
  "BG1" = 1452.782738, "BG2" = 426.3371158, "BG3" = 516.7650617,
  "BGS1" = 342.1892956, "BGS2" = 556.3127428, "BGS3" = 448.493284,
  "G1" = 93.29724454,  "G2" = 621.400516,  "G3" = 959.9439042,
  "GS1" = 381.1992513, "GS2" = 531.9358794, "GS3" = 1059.06864
)

# -----------------------------
# Robust date parser (no lubridate dependency)
# -----------------------------
.parse_dt <- function(x) {
  # Excel numeric serials -> seconds since 1899-12-30
  if (is.numeric(x)) return(as.POSIXct(x * 86400, origin = "1899-12-30", tz = "Europe/Paris"))
  # Try common text formats
  fmts <- c("%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S",
            "%m/%d/%Y %H:%M", "%d/%m/%Y %H:%M",
            "%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y")
  for (f in fmts) {
    dt <- suppressWarnings(as.POSIXct(x, format = f, tz = "Europe/Paris"))
    if (any(!is.na(dt))) return(dt)
  }
  # last resort
  suppressWarnings(as.POSIXct(x, tz = "Europe/Paris"))
}

# =============================
#       Load & Prepare Data
# =============================
cbd <- read.csv(in_csv, stringsAsFactors = FALSE)

# ---- detect and parse a date/time column into DateTime ----
date_candidates <- c("Date", "Datetime", "TIMESTAMP", "Date_Reco", "Date_GPP", "Date_AVAL")
date_col <- intersect(date_candidates, names(cbd))
if (length(date_col)) {
  cbd$DateTime <- .parse_dt(cbd[[date_col[1]]])
} else {
  warning("No obvious date column found; DateTime will be NA. Add one of: ",
          paste(date_candidates, collapse = ", "))
  cbd$DateTime <- NA
}

# Ensure expected columns are numeric
num_cols <- intersect(c("PAR","GPP","GPPmax_VI_pred","GPPmax_LAI_pred"), names(cbd))
cbd[num_cols] <- lapply(cbd[num_cols], function(x) suppressWarnings(as.numeric(x)))

cbd <- cbd %>%
  mutate(
    Replicates = as.character(Replicates),
    PlotID = factor(paste0(Plot, Replicates), levels = meso_order, ordered = TRUE),
    Veg    = factor(gsub("[0-9]+$", "", as.character(Plot)), levels = veg_levels, ordered = TRUE),
    K_value = as.numeric(k_values[as.character(PlotID)])
  )

# Modeled GPP using rectangular hyperbola with replicate-specific K
cbd <- cbd %>%
  mutate(
    GPP_VI_modeled  = ifelse(is.finite(GPPmax_VI_pred)  & is.finite(PAR) & is.finite(K_value),
                             (GPPmax_VI_pred  * PAR) / (K_value + PAR), NA_real_),
    GPP_LAI_modeled = ifelse(is.finite(GPPmax_LAI_pred) & is.finite(PAR) & is.finite(K_value),
                             (GPPmax_LAI_pred * PAR) / (K_value + PAR), NA_real_)
  )

# =============================
#    Long format + Metrics
# =============================
# Include DateTime in the long data so we can save it later
pieces_long <- list(
  cbd %>% transmute(DateTime, PlotID, Veg, Observed = GPP, Modeled = GPP_LAI_modeled, Model = "LAI")
)
if (USE_VI) {
  pieces_long <- c(pieces_long,
                   list(cbd %>% transmute(DateTime, PlotID, Veg, Observed = GPP, Modeled = GPP_VI_modeled,  Model = "VI")))
}

long_all <- bind_rows(pieces_long) %>%
  mutate(Model = factor(Model, levels = if (USE_VI) c("VI","LAI") else "LAI"))

# Helper: adjusted R^2 and p-value (overall F-test) safely
safe_adj_r2_p <- function(obs, pred){
  df <- data.frame(obs = obs, pred = pred) %>% tidyr::drop_na()
  if (nrow(df) < 3) return(list(r2 = NA_real_, p = NA_real_))
  fit <- tryCatch(lm(obs ~ pred, data = df), error = function(e) NULL)
  if (is.null(fit)) return(list(r2 = NA_real_, p = NA_real_))
  s <- summary(fit)
  r2 <- as.numeric(s$adj.r.squared)
  p  <- tryCatch({
    fs <- s$fstatistic
    if (is.null(fs)) NA_real_ else pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  }, error = function(e) NA_real_)
  list(r2 = r2, p = p)
}

# ===== Vectorized helpers for plot annotations =====
format_p_plotmath <- function(p){
  out <- ifelse(
    is.na(p), "italic(p)==NA",
    ifelse(p < 0.001, "italic(p)<0.001",
           paste0("italic(p)==", formatC(p, format = "f", digits = 3)))
  )
  as.character(out)
}
label_r2p <- function(model, r2, p) {
  model <- as.character(model)
  r2txt <- ifelse(is.na(r2), "NA", sprintf('%.2f', r2))
  ptxt  <- format_p_plotmath(p)
  ifelse(model == "VI",
         paste0("VI:~Adj~R^2==", r2txt, "~", ptxt),
         paste0("Adj~R^2==", r2txt, "~", ptxt))
}

metrics_repl <- long_all %>%
  group_by(Model, PlotID) %>%
  summarize(
    n     = sum(is.finite(Observed) & is.finite(Modeled)),
    stats = list(safe_adj_r2_p(Observed, Modeled)),
    .groups = "drop"
  ) %>%
  mutate(
    adjR2 = vapply(stats, function(z) z$r2, numeric(1)),
    Pval  = vapply(stats, function(z) z$p,  numeric(1))
  ) %>%
  select(-stats)

metrics_veg <- long_all %>%
  group_by(Model, Veg) %>%
  summarize(
    n     = sum(is.finite(Observed) & is.finite(Modeled)),
    stats = list(safe_adj_r2_p(Observed, Modeled)),
    .groups = "drop"
  ) %>%
  mutate(
    adjR2 = vapply(stats, function(z) z$r2, numeric(1)),
    Pval  = vapply(stats, function(z) z$p,  numeric(1))
  ) %>%
  select(-stats)

# Shared axis limits for 1:1 line and label placement
all_vals <- c(long_all$Observed, long_all$Modeled)
lim_min  <- floor(min(all_vals, na.rm = TRUE))
lim_max  <- ceiling(max(all_vals, na.rm = TRUE))

# =============================
# Helper: build horizontal legends
# =============================
.build_horizontal_legend <- function(){
  df <- data.frame(Model = factor(if (USE_VI) c("VI","LAI") else "LAI",
                                  levels = c("VI","LAI")),
                   x = 1, y = 1)
  p_leg <- ggplot(df, aes(x, y, fill = Model)) +
    geom_point(shape = 21, size = 4.5, stroke = 0.6, color = "black") +
    scale_fill_manual(values = c(VI = "black", LAI = "grey85"), name = NULL,
                      breaks = if (USE_VI) c("VI","LAI") else "LAI") +
    guides(fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(shape = 21, size = 4, color = "black", stroke = 0.6)
    )) +
    theme_void() +
    theme(
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.box.margin = margin(0,0,0,0),
      legend.margin     = margin(0,0,0,0),
      legend.spacing.x  = unit(12, "pt"),
      legend.text       = element_text(margin = margin(r = 10))
    )
  cowplot::get_legend(p_leg)
}
.build_horizontal_legend_no_lai <- function(){
  df <- data.frame(
    Model = factor(if (USE_VI) c("VI","LAI") else "LAI", levels = c("VI","LAI")),
    x = 1, y = 1
  )
  p_leg <- ggplot(df, aes(x, y, fill = Model)) +
    geom_point(shape = 21, size = 4.5, stroke = 0.6, color = "black") +
    scale_fill_manual(
      values  = c(VI = "black", LAI = "grey85"),
      breaks  = if (USE_VI) c("VI","LAI") else "LAI",
      labels  = if (USE_VI) c("VI", "") else c(""),  # remove "LAI" label
      name    = NULL
    ) +
    guides(fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(shape = 21, size = 4, color = "black", stroke = 0.6)
    )) +
    theme_void() +
    theme(
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.box.margin = margin(0,0,0,0),
      legend.margin     = margin(0,0,0,0),
      legend.spacing.x  = unit(12, "pt"),
      legend.text       = element_text(margin = margin(r = 10))
    )
  cowplot::get_legend(p_leg)
}

# =============================
# Plotting functions (unchanged)
# =============================
plot_model_combined_veg <- function(df_long, df_lab, file_out){
  rng   <- lim_max - lim_min
  xpad  <- 0.04 * rng
  ypad  <- 0.06 * rng
  gap   <- 0.08 * rng
  
  lab_df <- df_lab %>%
    mutate(
      adj_text = label_r2p(as.character(Model), adjR2, Pval),
      x = lim_min + xpad,
      y = lim_max - ypad - ifelse(Model == "LAI", gap, 0),
      hjust_lab = 0
    )
  
  title_txt <- if (USE_VI)
    "Observed vs Modeled GPP — VI (black) & LAI (light grey), pooled by vegetation"
  else
    "Observed vs Modeled GPP — LAI (light grey), pooled by vegetation"
  
  base_plot <- ggplot(df_long, aes(x = Observed, y = Modeled)) +
    geom_point(aes(fill = Model),
               shape = 21, color = "black", stroke = 0.6,
               alpha = 0.9, size = 2.1, na.rm = TRUE, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ Veg, ncol = 3, drop = FALSE) +
    geom_label(
      data = lab_df,
      aes(x = x, y = y, label = adj_text, hjust = hjust_lab),
      inherit.aes = FALSE, vjust = 1,
      color = "black", fill = "white", alpha = 0.9,
      label.size = 0.35, label.r = grid::unit(0.12, "lines"),
      size = 4, parse = TRUE
    ) +
    scale_x_continuous(limits = c(lim_min, lim_max),
                       expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(limits = c(lim_min, lim_max),
                       expand = expansion(mult = c(0.14, 0.14))) +
    scale_fill_manual(values = c(VI = "black", LAI = "grey85"),
                      breaks = if (USE_VI) c("VI","LAI") else "LAI",
                      guide = "none") +
    labs(
      title = title_txt,
      x = expression("Observed GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
      y = expression("Modeled GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")")
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 6)),
      legend.position = "none",
      plot.margin = margin(t = 6, r = 10, b = 10, l = 10)
    )
  
  leg <- .build_horizontal_legend_no_lai()
  p_final <- cowplot::plot_grid(leg, base_plot, ncol = 1, rel_heights = c(0.10, 1))
  ggsave(file_out, p_final, width = 12, height = 8, dpi = 300)
  message("Saved: ", file_out)
}

plot_model_combined_repl <- function(df_long, df_lab, file_out){
  rng   <- lim_max - lim_min
  xpad  <- 0.04 * rng
  ypad  <- 0.06 * rng
  gap   <- 0.08 * rng
  
  lab_df <- df_lab %>%
    mutate(
      adj_text = label_r2p(as.character(Model), adjR2, Pval),
      x = lim_min + xpad,
      y = lim_max - ypad - ifelse(Model == "LAI", gap, 0),
      hjust_lab = 0
    )
  
  title_txt <- if (USE_VI)
    "Observed vs Modeled GPP — VI (black) & LAI (light grey) per mesocosm"
  else
    "Observed vs Modeled GPP — LAI (light grey) per mesocosm"
  
  base_plot <- ggplot(df_long, aes(x = Observed, y = Modeled)) +
    geom_point(aes(fill = Model),
               shape = 21, color = "black", stroke = 0.6,
               alpha = 0.9, size = 1.9, na.rm = TRUE, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ PlotID, ncol = 3, drop = FALSE) +
    geom_label(
      data = lab_df,
      aes(x = x, y = y, label = adj_text, hjust = hjust_lab),
      inherit.aes = FALSE, vjust = 1,
      color = "black", fill = "white", alpha = 0.9,
      label.size = 0.35, label.r = grid::unit(0.12, "lines"),
      size = 3.6, parse = TRUE
    ) +
    scale_x_continuous(limits = c(lim_min, lim_max),
                       expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(limits = c(lim_min, lim_max),
                       expand = expansion(mult = c(0.14, 0.14))) +
    scale_fill_manual(values = c(VI = "black", LAI = "grey85"),
                      breaks = if (USE_VI) c("VI","LAI") else "LAI",
                      guide = "none") +
    labs(
      title = title_txt,
      x = expression("Observed GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
      y = expression("Modeled GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 6)),
      legend.position = "none",
      plot.margin = margin(t = 6, r = 10, b = 10, l = 10)
    )
  
  leg <- .build_horizontal_legend()
  p_final <- cowplot::plot_grid(leg, base_plot, ncol = 1, rel_heights = c(0.08, 1))
  ggsave(file_out, p_final, width = 10, height = 16, dpi = 300)
  message("Saved: ", file_out)
}

# =============================
#               RUN
# =============================
plot_model_combined_repl(
  df_long = long_all,
  df_lab  = metrics_repl,
  file_out = out_rep
)

plot_model_combined_veg(
  df_long = long_all,
  df_lab  = metrics_veg,
  file_out = out_veg
)

message("✅ Both figures written:\n- ", out_rep, "\n- ", out_veg)

# ==========================================================
#   ADDITIONAL: Time-series plots (Observed vs Modeled)
#   1) Replicate-level (3×5 facets)
#   2) Vegetation-mean ± SE (3×2 facets)
# ==========================================================

# ---- 1. Replicate-level long table ----
ts_long_repl <- cbd %>%
  select(DateTime, PlotID, Veg, GPP, GPP_LAI_modeled) %>%
  pivot_longer(cols = c(GPP, GPP_LAI_modeled),
               names_to = "Type", values_to = "Value") %>%
  mutate(Type = recode(Type,
                       GPP = "Observed",
                       GPP_LAI_modeled = "Modeled"),
         Type = factor(Type, levels = c("Observed", "Modeled")))

# ---- Plot 1: per replicate (3×5) ----
# ---- Plot 1: per replicate (3×5), POINTS ONLY ----
p_ts_repl <- ggplot(ts_long_repl,
                    aes(x = DateTime, y = Value, color = Type)) +
  geom_point(size = 1.2, alpha = 0.9, na.rm = TRUE) +
  facet_wrap(~ PlotID, ncol = 3, scales = "fixed") +
  scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
  labs(
    title = "Observed (black) vs Modeled (red) GPP — Replicates",
    x = "Date",
    y = expression("GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave("TimeSeries_GPP_byReplicate_LAI_3x5.png",
       p_ts_repl, width = 12, height = 15, dpi = 300)

# ---- 2. Compute vegetation means ± SE ----
se_fun <- function(x) {
  n <- sum(is.finite(x))
  if (n <= 1) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / sqrt(n)
}

veg_agg <- cbd %>%
  group_by(Veg, DateTime) %>%
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
  mutate(Type = recode(Type,
                       Obs = "Observed",
                       Mod = "Modeled"),
         Type = factor(Type, levels = c("Observed", "Modeled")))

# ---- Plot 2: mean ± SE per vegetation (3×2) ----
# ---- Plot 2: mean ± SE per vegetation (3×2), POINTS + ERRORBARS ----
p_ts_veg <- ggplot(veg_agg,
                   aes(x = DateTime, y = mean, color = Type)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0, alpha = 0.7, na.rm = TRUE) +
  geom_point(size = 1.3, alpha = 0.9, na.rm = TRUE) +
  facet_wrap(~ Veg, ncol = 3, scales = "fixed") +
  scale_color_manual(values = c("Observed" = "black", "Modeled" = "red")) +
  labs(
    title = "Observed vs Modeled GPP (mean ± SE) — Vegetation Groups",
    x = "Date",
    y = expression("Mean GPP ("*mu*"mol CO"[2]*" m"^-2*" s"^-1*")"),
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave("TimeSeries_GPP_byVeg_LAI_3x2.png",
       p_ts_veg, width = 12, height = 8, dpi = 300)

# =============================
# Save GPP Observed vs Modeled (Evaluation sheet) WITH DATE
# =============================
save_eval_workbook <- function(df_long, out_file){
  # Keep only LAI modeled version (since USE_VI = FALSE)
  df_eval <- df_long %>%
    filter(Model == "LAI") %>%
    transmute(
      Veg,
      PlotID,
      Date = ifelse(is.na(DateTime), NA_character_,
                    format(as.POSIXct(DateTime, origin = "1970-01-01", tz = "Europe/Paris"),
                           "%Y-%m-%d %H:%M:%S")),
      GPP_obs = Observed,
      GPP_mod = Modeled,
      Residual = Observed - Modeled
    ) %>%
    arrange(Veg, PlotID, Date)
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Evaluation")
  openxlsx::writeData(wb, "Evaluation", df_eval)
  
  # Optional placeholders to mirror other files
  openxlsx::addWorksheet(wb, "Calibration")
  openxlsx::writeData(wb, "Calibration",
                      data.frame(note = "Calibration results not included here"))
  openxlsx::addWorksheet(wb, "Combined")
  openxlsx::writeData(wb, "Combined",
                      data.frame(note = "Combined outputs not included here"))
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  message("✔ Wrote Evaluation workbook (with Date): ", out_file)
}

# =============================
# Call after plots
# =============================
eval_outfile <- "GPP_Observed_vs_Modeled_All_LAI.xlsx"
save_eval_workbook(long_all, eval_outfile)

# ===========================================================
#   Compute adj R², RMSE, and model p-value for each replicate
#   from GPP_Observed_vs_Modeled_All_LAI (with correct columns)
# ===========================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(hydroGOF)
library(readxl)
library(openxlsx)

# === 1. Load data ===
file_path <- "GPP_Observed_vs_Modeled_All_LAI.xlsx"

if (grepl("\\.csv$", file_path)) {
  gpp_df <- read.csv(file_path)
} else {
  gpp_df <- read_excel(file_path)
}

# Check columns
print(names(gpp_df))

# === 2. Assign correct names ===
rep_col <- "PlotID"
obs_col <- "GPP_obs"
mod_col <- "GPP_mod"

# === 3. Compute metrics per replicate ===
calc_metrics <- function(df, obs_col, mod_col) {
  y_obs <- df[[obs_col]]
  y_mod <- df[[mod_col]]
  keep <- is.finite(y_obs) & is.finite(y_mod)
  if (sum(keep) < 3) return(c(R2 = NA, RMSE = NA, pval_model = NA))
  
  lm_fit <- lm(y_obs[keep] ~ y_mod[keep])
  adjR2 <- summary(lm_fit)$adj.r.squared
  RMSE  <- hydroGOF::rmse(y_obs[keep], y_mod[keep])
  
  # p-value for full model (ANOVA F-test)
  aov_sum <- anova(lm_fit)
  F_val <- aov_sum$`F value`[1]
  df1 <- aov_sum$Df[1]; df2 <- aov_sum$Df[2]
  pval_model <- 1 - pf(F_val, df1, df2)
  
  c(R2 = adjR2, RMSE = RMSE, pval_model = pval_model)
}

metrics_df <- gpp_df %>%
  group_by(!!sym(rep_col)) %>%
  summarise(
    R2 = calc_metrics(cur_data(), obs_col, mod_col)[["R2"]],
    RMSE = calc_metrics(cur_data(), obs_col, mod_col)[["RMSE"]],
    pval_model = calc_metrics(cur_data(), obs_col, mod_col)[["pval_model"]],
    .groups = "drop"
  ) %>%
  mutate(SigFlag = ifelse(pval_model > 0.05, "ns", "sig"))

# === 4. Save metrics table ===
write.xlsx(metrics_df, "GPP_Model_Performance_byReplicate.xlsx", rowNames = FALSE)
cat("✅ Saved: GPP_Model_Performance_byReplicate.xlsx\n")

# === 5. Plot RMSE vs Adjusted R² ===
ggplot(metrics_df, aes(x = R2, y = RMSE, label = !!sym(rep_col), color = SigFlag)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 4, max.overlaps = 100, box.padding = 0.3, segment.color = "gray60") +
  scale_color_manual(values = c("sig" = "black", "ns" = "red"),
                     name = "Model significance (p > 0.05)") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 14)) +
  theme_bw(base_size = 18) +
  labs(
    title = "Calibration RMSE vs Adjusted R² per Replicate\n(red = non-significant fits)",
    x = expression(Adjusted~R^2),
    y = "RMSE"
  ) +
  theme(strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 14),
        legend.position = "top")

ggsave("Calibration_RMSE_vs_R2_perReplicate_from_GPP_All_LAI.png",
       width = 12, height = 8, dpi = 300)

cat("✅ Saved: Calibration_RMSE_vs_R2_perReplicate_from_GPP_All_LAI.png (red = non-significant)\n")
