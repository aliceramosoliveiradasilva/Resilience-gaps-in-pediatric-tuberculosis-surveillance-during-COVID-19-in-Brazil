# ============================================================
# Multi-panel ITS for TB Cases by Ethnicity (A–F)
# - Shared legend; pandemic shaded (>=2020)
# - Model fit = dashed (thicker); Counterfactual = solid (thinner)
# - Exports PDF, PNG, SVG
# ============================================================

rm(list = ls())
if (!is.null(dev.list())) try(dev.off(), silent = TRUE)

# Packages
pkgs <- c("tidyverse","janitor","broom","lmtest","sandwich","nlme","patchwork","stringi")
new <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(pkgs, library, character.only = TRUE))

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b

# ---------- paths ----------
setwd("/home/ramos/Documentos/R/Subnotificacao TB/")
out_dir <- "saidas_its"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ---------- read & standardize ----------
dados <- read.csv("Taxaa_longo.csv", sep = ";", check.names = FALSE) |>
  janitor::clean_names()

stopifnot(all(c("ano","taxa_100k","etnia","indicador") %in% names(dados)))

map_eth_pt2en <- function(x){
  x_low <- stringr::str_to_lower(trimws(x))
  dplyr::case_when(
    x_low %in% c("branca","branco")                  ~ "White",
    x_low %in% c("preta","preto")                    ~ "Black",
    x_low %in% c("parda","pardo")                    ~ "Brown (Pardo)",
    x_low %in% c("amarelo","amarela")                ~ "Asian",
    x_low %in% c("indigena","indígena","indiginea")  ~ "Indigenous",
    x_low %in% c("total")                            ~ "Total",
    TRUE                                             ~ stringr::str_to_title(x)
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
    etnia = map_eth_pt2en(etnia),
    indicador = map_ind_pt2en(indicador)
  ) |>
  filter(is.finite(ano), is.finite(taxa_100k))

# ---------- keep only Cases ----------
dados_cases <- dados |> filter(indicador == "Cases")

# ---------- ITS helpers ----------
make_its_vars <- function(df, break_year = 2020){
  df <- df |> arrange(ano)
  df |>
    mutate(
      t = row_number(),
      post = if_else(ano >= break_year, 1L, 0L),
      t_post = if_else(post == 1L, t - min(t[post == 1L], na.rm = TRUE) + 1L, 0L)
    )
}

fit_its_models <- function(df, y = "taxa_100k", log_scale = TRUE){
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
  
  # GLS AR(1)
  gls_fit <- try(
    nlme::gls(f, data = d, method = "ML",
              correlation = nlme::corARMA(form = ~ t, p = 1)),
    silent = TRUE
  )
  
  # escolha
  lb_ols <- Box.test(residuals(ols), type = "Ljung-Box", lag = min(10, Tn-1))
  dw_ols <- lmtest::dwtest(ols)
  has_ac <- (lb_ols$p.value %||% 1) < 0.05 || (dw_ols$p.value %||% 1) < 0.05
  
  if (!inherits(gls_fit, "try-error")) {
    aic_ols <- AIC(ols)
    aic_gls <- AIC(gls_fit)
    pick <- if (has_ac) "GLS" else if (is.finite(aic_gls) && aic_gls < aic_ols) "GLS" else "OLS"
  } else {
    pick <- "OLS"
  }
  
  list(
    data = d,
    ols = list(fit = ols, vcov = vcov_nw),
    gls = list(fit = if (inherits(gls_fit, "try-error")) NULL else gls_fit),
    chosen = pick,
    log_scale = log_scale
  )
}

make_fitted_cf <- function(res, break_year = 2020){
  d <- res$data
  cf <- if (res$chosen == "GLS" && !is.null(res$gls$fit)) coef(res$gls$fit) else coef(res$ols$fit)
  with(d, {
    y_hat <- cf[1] + cf["t"]*t + cf["post"]*post + cf["t_post"]*t_post
    y_cf  <- cf[1] + cf["t"]*t
    tibble(ano, y_hat = y_hat, y_cf = y_cf)
  }) |>
    mutate(
      y_hat = if (res$log_scale) exp(y_hat) else y_hat,
      y_cf  = if (res$log_scale) exp(y_cf)  else y_cf
    )
}

# ---------- which ethnicities and facet labels (A–F) ----------
facet_order <- c("Total","Asian","Black","Brown (Pardo)","White","Indigenous")
facet_label <- c(
  "Total"         = "A - Total",
  "Asian"         = "B - Asian",
  "Black"         = "C - Black",
  "Brown (Pardo)" = "D - Brown (Pardo)",
  "White"         = "E - White",
  "Indigenous"    = "F - Indigenous"
)

avail <- dados_cases |> count(etnia, name = "n") |> filter(etnia %in% facet_order, n >= 4) |> pull(etnia)
eth_to_plot <- facet_order[facet_order %in% avail]

# ---------- build plotting data ----------
break_year <- 2020
all_plot <- list()

for (eth in eth_to_plot) {
  df <- dados_cases |>
    filter(etnia == eth) |>
    arrange(ano) |>
    make_its_vars(break_year = break_year)
  
  if (nrow(df) < 4) next
  
  fit <- fit_its_models(df, y = "taxa_100k", log_scale = TRUE)
  fc  <- make_fitted_cf(fit, break_year = break_year)
  
  tmp <- df |>
    select(ano, taxa_100k) |>
    mutate(Series = "Observed") |>
    rename(value = taxa_100k)
  
  model_line <- fc |>
    transmute(ano, value = y_hat, Series = "Model (fitted)")
  
  cf_line <- fc |>
    transmute(ano, value = y_cf, Series = "Counterfactual")
  
  all_plot[[length(all_plot)+1]] <- bind_rows(tmp, model_line, cf_line) |>
    mutate(Ethnicity = eth,
           Facet = facet_label[eth])
}

plot_df <- bind_rows(all_plot) |>
  mutate(
    Facet = factor(Facet, levels = facet_label[eth_to_plot]),
    Series = factor(Series, levels = c("Observed","Model (fitted)","Counterfactual"))
  )

# ---------- plotting ----------
cols  <- c("Observed" = "#222222", "Model (fitted)" = "#1f77b4", "Counterfactual" = "#d62728")
lines <- c("Observed" = "solid",   "Model (fitted)" = "longdash", "Counterfactual" = "solid")
sizes <- c("Observed" = 0.6,       "Model (fitted)" = 0.9,        "Counterfactual" = 0.5)

x_min <- min(dados_cases$ano, na.rm = TRUE)
x_max <- max(dados_cases$ano, na.rm = TRUE)

p <- ggplot() +
  annotate("rect",
           xmin = break_year - 0.5, xmax = x_max + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "#444444", alpha = 0.20) +
  geom_point(data = plot_df |> filter(Series == "Observed"),
             aes(x = ano, y = value),
             color = cols["Observed"], size = 1.8, alpha = 0.9, show.legend = FALSE) +
  geom_line(data = plot_df |> filter(Series == "Observed"),
            aes(x = ano, y = value, color = Series, linetype = Series, size = Series)) +
  geom_line(data = plot_df |> filter(Series == "Model (fitted)"),
            aes(x = ano, y = value, color = Series, linetype = Series, size = Series)) +
  geom_line(data = plot_df |> filter(Series == "Counterfactual"),
            aes(x = ano, y = value, color = Series, linetype = Series, size = Series)) +
  facet_wrap(~ Facet, ncol = 3, scales = "free_y") +
  scale_color_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = lines, name = NULL) +
  scale_size_manual(values = sizes, guide = guide_legend(override.aes = list(alpha = 1)), name = NULL) +
  scale_x_continuous(breaks = seq(x_min, x_max, by = 2), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Year", y = "Rate per 100,000") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank()
  )

# ---------- export ----------
stub <- file.path(out_dir, "FIG_Cases_by_Ethnicity_A_to_F")
ggsave(paste0(stub, ".png"), p, width = 12, height = 8, dpi = 300)
ggsave(paste0(stub, ".pdf"), p, width = 12, height = 8, device = "pdf")

if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
ggsave(paste0(stub, ".svg"), p, width = 12, height = 8, device = svglite::svglite)

message("Figure exported:\n",
        " - ", paste0(stub, ".png"), "\n",
        " - ", paste0(stub, ".pdf"), "\n",
        " - ", paste0(stub, ".svg"), "\n")


