# ============================================================
# ITS by Ethnicity x Indicator (break in 2020)
# - OLS + Newey-West and GLS AR(1)
# - Auto-selection of best model per series
# - Ethnicity & indicator names in English
# - Pandemic (>=2020) shaded area
# - Year ticks shown for every year
# - Export PNG and PDF
# ============================================================

# -------------------------
# 0) Setup
# -------------------------
rm(list = ls())
if (!is.null(dev.list)) try(dev.off(), silent = TRUE)

pkgs <- c("tidyverse","janitor","broom","lmtest","sandwich","nlme","zoo","patchwork","stringi")
new <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(pkgs, library, character.only = TRUE))

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b

# -------------------------
# 1) Read & standardize
# -------------------------
setwd("/home/ramos/Documentos/R/Subnotificacao TB/")

dados <- read.csv("Taxaa_longo.csv", sep = ";", check.names = FALSE) |>
  janitor::clean_names()

stopifnot(all(c("ano","taxa_100k","etnia","indicador") %in% names(dados)))

# Map Portuguese -> English for Ethnicity and Indicator
map_eth_pt2en <- function(x){
  x_low <- stringr::str_to_lower(trimws(x))
  dplyr::case_when(
    x_low %in% c("branca","branco")               ~ "White",
    x_low %in% c("preta","preto")                 ~ "Black",
    x_low %in% c("parda","pardo")                 ~ "Brown (Pardo)",
    x_low %in% c("amarelo","amarela")             ~ "Asian",
    x_low %in% c("indigena","indígena","indiginea") ~ "Indigenous",
    x_low %in% c("total")                         ~ "Total",
    TRUE                                          ~ stringr::str_to_title(x)
  )
}

map_ind_pt2en <- function(x){
  x_low <- stringr::str_to_lower(trimws(x))
  dplyr::case_when(
    stringr::str_detect(x_low, "^cas")  ~ "Cases",
    stringr::str_detect(x_low, "^obit") ~ "Deaths",
    stringr::str_detect(x_low, "^inte") ~ "Hospitalizations",
    TRUE                                ~ stringr::str_to_title(x)
  )
}

dados <- dados |>
  mutate(
    etnia       = map_eth_pt2en(etnia),
    indicador   = map_ind_pt2en(indicador)
  ) |>
  arrange(indicador, etnia, ano) |>
  filter(is.finite(ano), is.finite(taxa_100k))

# -------------------------
# 2) Helpers (ITS & models)
# -------------------------

# ITS variables
make_its_vars <- function(df, break_year = 2020){
  df <- df |> arrange(ano)
  df |>
    mutate(
      t = row_number(),
      post = if_else(ano >= break_year, 1L, 0L),
      t_post = if_else(post == 1L, t - min(t[post == 1L], na.rm = TRUE) + 1L, 0L)
    )
}

# Tidy/Glance for nlme::gls
tidy_gls_nlme <- function(model, conf.int = TRUE, conf.level = 0.95){
  sm <- summary(model)
  tt <- as.data.frame(sm$tTable)
  tt$term <- rownames(sm$tTable)
  names(tt) <- sub("^Value$", "estimate", names(tt))
  names(tt) <- sub("^Std\\.Error$", "std.error", names(tt))
  names(tt) <- sub("^t-value$", "statistic", names(tt))
  names(tt) <- sub("^p-value$", "p.value", names(tt))
  tt <- tt[, c("term","estimate","std.error","statistic","p.value")]
  
  if (conf.int) {
    ci <- try(nlme::intervals(model, which = "coef", level = conf.level), silent = TRUE)
    if (!inherits(ci, "try-error")) {
      ci_tab <- as.data.frame(ci$coef)
      ci_tab$term <- rownames(ci$coef)
      names(ci_tab) <- c("conf.low","estimate_ci","conf.high","term")
      tt <- dplyr::left_join(tt, ci_tab[, c("term","conf.low","conf.high")], by = "term")
    } else {
      z <- qnorm(0.5 + conf.level/2)
      tt$conf.low  <- tt$estimate - z*tt$std.error
      tt$conf.high <- tt$estimate + z*tt$std.error
    }
  }
  tibble::as_tibble(tt)
}

glance_gls_nlme <- function(model){
  tibble::tibble(
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    nobs = nobs(model)
  )
}

# Fit OLS+NW and GLS AR(1); optional log scale
fit_its_models <- function(df, y = "taxa_100k", log_scale = FALSE){
  d <- df
  if (log_scale) {
    min_pos <- suppressWarnings(min(d[[y]][d[[y]] > 0], na.rm = TRUE))
    eps <- ifelse(is.finite(min_pos), min_pos*0.01, 1e-6)
    d <- d |> mutate(y_var = log(pmax(.data[[y]], eps)))
  } else {
    d <- d |> mutate(y_var = .data[[y]])
  }
  
  f <- as.formula("y_var ~ t + post + t_post")
  
  # OLS + Newey-West
  ols <- lm(f, data = d)
  Tn <- nrow(d)
  lag_nw <- max(0, floor(Tn^(1/4)))
  vcov_nw <- sandwich::NeweyWest(ols, lag = lag_nw, prewhite = FALSE)
  tidy_ols <- broom::tidy(ols, conf.int = TRUE, conf.level = 0.95, vcov = vcov_nw)
  glance_ols <- broom::glance(ols)
  rmse_ols <- sqrt(mean(residuals(ols)^2))
  dw_ols <- lmtest::dwtest(ols)
  lb_ols <- Box.test(residuals(ols), type = "Ljung-Box", lag = min(10, Tn-1))
  
  # GLS AR(1)
  gls_fit <- try(
    nlme::gls(f, data = d, method = "ML", correlation = nlme::corARMA(form = ~ t, p = 1)),
    silent = TRUE
  )
  
  if (!inherits(gls_fit, "try-error")) {
    tidy_gls <- tidy_gls_nlme(gls_fit, conf.int = TRUE, conf.level = 0.95)
    glance_gls <- glance_gls_nlme(gls_fit)
    rmse_gls <- sqrt(mean(residuals(gls_fit, type = "response")^2))
    res_gls <- residuals(gls_fit, type = "response")
    lb_gls <- Box.test(res_gls, type = "Ljung-Box", lag = min(10, Tn-1))
  } else {
    tidy_gls <- tibble::tibble(term = c("(Intercept)","t","post","t_post"),
                               estimate = NA_real_, std.error = NA_real_,
                               statistic = NA_real_, p.value = NA_real_,
                               conf.low = NA_real_, conf.high = NA_real_)
    glance_gls <- tibble::tibble(AIC = NA_real_, BIC = NA_real_, logLik = NA_real_, nobs = nrow(d))
    rmse_gls <- NA_real_
    lb_gls <- NA
  }
  
  list(
    data = d,
    ols = list(
      fit = ols, tidy = tidy_ols, glance = glance_ols, rmse = rmse_ols,
      dw = dw_ols, lb = lb_ols, lag_nw = lag_nw, log_scale = log_scale
    ),
    gls = list(
      fit = if (inherits(gls_fit, "try-error")) NULL else gls_fit,
      tidy = tidy_gls, glance = glance_gls, rmse = rmse_gls,
      lb = lb_gls, log_scale = log_scale
    )
  )
}

# Choose best model using autocorrelation + AIC
choose_best_model <- function(res_lin_or_log){
  ols_lb_p <- res_lin_or_log$ols$lb$p.value %||% 1
  ols_dw_p <- res_lin_or_log$ols$dw$p.value %||% 1
  has_ac_ols <- (ols_lb_p < 0.05) || (ols_dw_p < 0.05)
  gls_avail <- !is.null(res_lin_or_log$gls$fit)
  
  if (has_ac_ols && gls_avail) return("GLS")
  if (gls_avail) {
    aic_ols <- res_lin_or_log$ols$glance$AIC %||% Inf
    aic_gls <- res_lin_or_log$gls$glance$AIC %||% Inf
    return(ifelse(aic_gls < aic_ols, "GLS", "OLS"))
  }
  "OLS"
}

# Fitted & counterfactual
make_fitted_cf <- function(res_list, break_year = 2020, use_model = c("auto","OLS","GLS")){
  use_model <- match.arg(use_model)
  d <- res_list$data
  
  pick <- function(){
    if (use_model == "OLS") return(coef(res_list$ols$fit))
    if (use_model == "GLS" && !is.null(res_list$gls$fit)) return(coef(res_list$gls$fit))
    best <- choose_best_model(res_list)
    if (best == "GLS") return(coef(res_list$gls$fit))
    coef(res_list$ols$fit)
  }
  
  cf <- pick()
  with(d, {
    y_hat <- cf[1] + cf["t"]*t + cf["post"]*post + cf["t_post"]*t_post
    y_cf  <- cf[1] + cf["t"]*t
    tibble::tibble(ano, t, post, t_post, y_hat = y_hat, y_cf = y_cf)
  })
}

# Effects table (post = level; t_post = slope)
effect_table <- function(res_list, ethnicity, indicator, scale_label, chosen){
  use_gls <- (chosen == "GLS") && !is.null(res_list$gls$fit)
  tt <- if (use_gls) res_list$gls$tidy else res_list$ols$tidy
  
  tt |>
    dplyr::filter(term %in% c("post","t_post")) |>
    dplyr::transmute(
      Ethnicity = ethnicity,
      Indicator = indicator,
      Scale = scale_label,
      Chosen_model = chosen,
      Effect = dplyr::recode(term,
                             post = "Change in level (post-2020)",
                             t_post = "Change in slope (post-2020)"),
      Estimate = estimate,
      CI95_low = conf.low,
      CI95_high = conf.high,
      p_value = p.value
    )
}

# Sanitize file names (ASCII only, underscores)
sanitize_stub <- function(text){
  x <- stringi::stri_trans_general(text, "Latin-ASCII")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+$", "", x)
  x
}

# Plot ITS (English, pandemic shaded, year ticks each year)
plot_its <- function(res_list, serie_label, break_year = 2020, original_y, log_scale, chosen){
  d <- res_list$data
  fits <- make_fitted_cf(res_list, break_year = break_year, use_model = "auto")
  
  # Back-transform when using log-scale
  if (log_scale) {
    d$y_plot <- exp(d$y_var)
    fits <- fits |> dplyr::mutate(y_hat = exp(y_hat), y_cf = exp(y_cf))
    subttl <- paste0("Log scale (~percent changes). Chosen model: ", chosen)
  } else {
    d$y_plot <- d$y_var
    subttl <- paste0("Original scale. Chosen model: ", chosen)
  }
  
  x_min <- min(d$ano, na.rm = TRUE)
  x_max <- max(d$ano, na.rm = TRUE)
  
  ggplot(d, aes(x = ano, y = y_plot)) +
    # Pandemic shading (>= 2020)
    annotate("rect",
             xmin = break_year - 0.5, xmax = x_max + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "grey70", alpha = 0.2) +
    geom_point() +
    geom_line() +
    geom_line(data = fits, aes(x = ano, y = y_cf), linetype = 2) +
    geom_line(data = fits, aes(x = ano, y = y_hat)) +
    geom_vline(xintercept = break_year, linetype = 3) +
    scale_x_continuous(breaks = seq(x_min, x_max, by = 1)) +
    labs(
      title = paste0("Interrupted Time Series — ", serie_label),
      subtitle = subttl,
      x = "Year", y = original_y,
      caption = "Shaded area: pandemic period (>= 2020). Dashed line: counterfactual (projection of pre-2020 trend)."
    ) +
    theme_minimal(base_size = 12)
}

# Save PNG and PDF
save_both <- function(file_stub, plot, width = 8, height = 5, dpi = 300){
  ggsave(paste0(file_stub, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(file_stub, ".pdf"), plot, width = width, height = height, device = "pdf")
}

# -------------------------
# 3) Batch run
# -------------------------
out_dir <- "saidas_its"
if (!dir.exists(out_dir)) dir.create(out_dir)

break_year <- 2020

series <- dados |>
  dplyr::count(etnia, indicador, name = "n") |>
  dplyr::filter(n >= 4)

effects_all <- list()
quality_all <- list()

for (i in seq_len(nrow(series))) {
  et <- series$etnia[i]
  ind <- series$indicador[i]
  
  df <- dados |>
    dplyr::filter(etnia == et, indicador == ind) |>
    dplyr::arrange(ano) |>
    make_its_vars(break_year = break_year)
  
  if (nrow(df) < 4) next
  
  # Both scales
  res_lin <- fit_its_models(df, y = "taxa_100k", log_scale = FALSE)
  res_log <- fit_its_models(df, y = "taxa_100k", log_scale = TRUE)
  
  # Auto-choose per scale
  chosen_lin <- choose_best_model(res_lin)
  chosen_log <- choose_best_model(res_log)
  
  # Effects (chosen model per scale)
  effects_all[[length(effects_all)+1]] <- effect_table(res_lin, et, ind, "Linear", chosen_lin)
  effects_all[[length(effects_all)+1]] <- effect_table(res_log, et, ind, "Log",    chosen_log)
  
  # Quality / diagnostics
  qa <- tibble::tibble(
    Ethnicity = et, Indicator = ind,
    Scale = c("Linear","Linear","Log","Log"),
    Model = c("OLS+NW","GLS AR(1)","OLS+NW","GLS AR(1)"),
    Chosen_model = c(chosen_lin, chosen_lin, chosen_log, chosen_log),
    RMSE = c(res_lin$ols$rmse, res_lin$gls$rmse, res_log$ols$rmse, res_log$gls$rmse),
    AIC  = c(res_lin$ols$glance$AIC %||% NA_real_,
             res_lin$gls$glance$AIC %||% NA_real_,
             res_log$ols$glance$AIC %||% NA_real_,
             res_log$gls$glance$AIC %||% NA_real_),
    BIC  = c(res_lin$ols$glance$BIC %||% NA_real_,
             res_lin$gls$glance$BIC %||% NA_real_,
             res_log$ols$glance$BIC %||% NA_real_,
             res_log$gls$glance$BIC %||% NA_real_),
    DW_p = c(res_lin$ols$dw$p.value %||% NA_real_, NA, res_log$ols$dw$p.value %||% NA_real_, NA),
    LB_p = c(res_lin$ols$lb$p.value %||% NA_real_,
             if (is.list(res_lin$gls$lb)) res_lin$gls$lb$p.value else res_lin$gls$lb %||% NA_real_,
             res_log$ols$lb$p.value %||% NA_real_,
             if (is.list(res_log$gls$lb)) res_log$gls$lb$p.value else res_log$gls$lb %||% NA_real_)
  )
  quality_all[[length(quality_all)+1]] <- qa
  
  # Labels and file stubs
  lbl <- paste(et, "—", ind)
  stub_lin <- file.path(out_dir, paste0("ITS_", sanitize_stub(et), "_", sanitize_stub(ind), "_linear"))
  stub_log <- file.path(out_dir, paste0("ITS_", sanitize_stub(et), "_", sanitize_stub(ind), "_log"))
  
  # Plots (English, shaded pandemic, annual ticks)
  p_lin <- plot_its(res_lin, lbl, break_year, original_y = "Rate per 100,000", log_scale = FALSE, chosen = chosen_lin)
  p_log <- plot_its(res_log, lbl, break_year, original_y = "Rate per 100,000 (original scale)", log_scale = TRUE, chosen = chosen_log)
  
  # Save PNG & PDF
  save_both(stub_lin, p_lin, width = 8, height = 5)
  save_both(stub_log, p_log, width = 8, height = 5)
}

# Consolidate and export tables
effects_tbl <- dplyr::bind_rows(effects_all)
quality_tbl <- dplyr::bind_rows(quality_all)

readr::write_csv(effects_tbl, file.path(out_dir, "ITS_effects_by_series.csv"))
readr::write_csv(quality_tbl, file.path(out_dir, "ITS_model_quality.csv"))

# -------------------------
# 4) Total panel (English + shaded)
# -------------------------
if ("Total" %in% unique(dados$etnia)) {
  df_total <- dados |>
    dplyr::filter(etnia == "Total") |>
    dplyr::group_by(indicador) |>
    dplyr::group_modify(~ make_its_vars(.x, break_year = break_year)) |>
    dplyr::ungroup()
  
  plots <- df_total |>
    dplyr::group_split(indicador) |>
    purrr::map(function(df_i){
      res <- fit_its_models(df_i, y = "taxa_100k", log_scale = FALSE)
      chosen <- choose_best_model(res)
      plot_its(res,
               serie_label = paste0("Total — ", unique(df_i$indicador)),
               break_year = break_year,
               original_y = "Rate per 100,000",
               log_scale = FALSE,
               chosen = chosen)
    })
  
  if (length(plots) > 0) {
    p_all <- patchwork::wrap_plots(plots, ncol = 1)
    save_both(file.path(out_dir, "ITS_Total_panel_linear"), p_all, width = 8, height = 5*length(plots))
  }
}

# -------------------------
# 5) Done
# -------------------------
message("\nDone!\n",
        "- Output directory: ", normalizePath(out_dir), "\n",
        "- Tables:\n   * ITS_effects_by_series.csv (chosen model per scale)\n",
        "   * ITS_model_quality.csv (RMSE, AIC/BIC, DW/LB, chosen model)\n",
        "- Figures: PNG and PDF per series (linear & log) + Total panel.\n",
        "- Ethnicity & indicators labeled in EN; pandemic (>=2020) shaded; year ticks every year.\n")


