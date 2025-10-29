# ============================================================
# ITS -> Tabela .DOCX (editável) com flextable + officer
# - Todos os valores numéricos formatados com 2 casas decimais
# ============================================================

# ---- Pacotes ----
pkgs <- c("dplyr","readr","flextable","officer","stringr","tidyr")
new <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Caminhos ----
dir_out  <- "saidas_its"
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
arq_csv  <- file.path(dir_out, "ITS_effects_by_series.csv")
arq_docx <- file.path(dir_out, "Table1_ITS_effects.docx")

# ---- Ler dados ----
effects_tbl <- readr::read_csv(arq_csv, show_col_types = FALSE)

# ---- Recodificar etnias (se vierem em PT) ----
recode_eth <- function(x){
  x2 <- as.character(x)
  dplyr::recode(x2,
                "Branca"    = "White",
                "Preta"     = "Black",
                "Parda"     = "Brown/Mixed",
                "Amarela"   = "Asian",
                "Indígena"  = "Indigenous",
                .default    = x2
  )
}

# ---- Filtrar escala Log e converter em % ----
tab_resumo <- effects_tbl %>%
  dplyr::filter(Scale == "Log") %>%
  dplyr::mutate(
    Ethnicity = recode_eth(Ethnicity),
    Effect = dplyr::recode(
      Effect,
      "Change in level (post-2020)" = "Immediate level change (2020)",
      "Change in slope (post-2020)" = "Trend change after 2020",
      .default = Effect
    ),
    Percent = (exp(Estimate) - 1) * 100,
    CI_low  = (exp(CI95_low)  - 1) * 100,
    CI_high = (exp(CI95_high) - 1) * 100,
    p_value = as.numeric(p_value)
  ) %>%
  dplyr::select(Ethnicity, Indicator, Effect, Percent, CI_low, CI_high, p_value, Chosen_model) %>%
  dplyr::arrange(Ethnicity, Indicator, Effect)

# ---- Coluna concatenada Estimate (95% CI) ----
tab_resumo <- tab_resumo %>%
  dplyr::mutate(
    Estimate_CI = sprintf("%.2f%% (%.2f%%; %.2f%%)", Percent, CI_low, CI_high)
  )

# ---- Asteriscos de significância ----
p_to_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05, "*", "")
                )
         )
  )
}
tab_resumo <- tab_resumo %>%
  dplyr::mutate(sig = p_to_stars(p_value))

# ---- Escolha formato da tabela: A (colunas separadas) ou B (Estimate + IC juntos) ----
mostrar_forma <- "A"  # use "A" ou "B"

if(mostrar_forma == "A"){
  df_show <- tab_resumo %>%
    dplyr::select(Ethnicity, Indicator, Effect, Percent, CI_low, CI_high,
                  p_value, sig, Chosen_model)
} else {
  df_show <- tab_resumo %>%
    dplyr::select(Ethnicity, Indicator, Effect, Estimate_CI, p_value, sig, Chosen_model)
}

ft <- flextable::flextable(df_show)

# ---- Rótulos e formatação ----
if(mostrar_forma == "A"){
  ft <- flextable::set_header_labels(
    ft,
    Ethnicity     = "Ethnicity",
    Indicator     = "Outcome",
    Effect        = "Effect",
    Percent       = "Estimate",
    CI_low        = "95% CI (low)",
    CI_high       = "95% CI (high)",
    p_value       = "p-value",
    sig           = "",
    Chosen_model  = "Model"
  )
  ft <- flextable::add_header_row(
    ft,
    values = c("", "", "", "", "95% Confidence Interval", "", "", ""),
    colwidths = c(1, 1, 1, 1, 2, 1, 1, 1)
  )
  # Format numérico: 2 casas + %
  ft <- flextable::colformat_num(
    ft, j = c("Percent","CI_low","CI_high"),
    digits = 2, suffix = "%"
  )
} else {
  ft <- flextable::set_header_labels(
    ft,
    Ethnicity     = "Ethnicity",
    Indicator     = "Outcome",
    Effect        = "Effect",
    Estimate_CI   = "Estimate (95% CI)",
    p_value       = "p-value",
    sig           = "",
    Chosen_model  = "Model"
  )
}

# p-value com 2 casas
ft <- flextable::colformat_num(ft, j = "p_value", digits = 2)

# ---- Aparência ----
ft <- flextable::merge_v(ft, j = "Ethnicity")
ft <- flextable::autofit(ft)
ft <- flextable::theme_booktabs(ft)

col_left  <- intersect(c("Indicator","Effect","Chosen_model","Estimate_CI"), names(df_show))
if(length(col_left)) ft <- flextable::align(ft, j = col_left, align = "left", part = "body")

col_center <- intersect(c("Percent","CI_low","CI_high","p_value","sig"), names(df_show))
if(length(col_center)) ft <- flextable::align(ft, j = col_center, align = "center", part = "body")

ft <- flextable::align(ft, part = "header", align = "center")
ft <- flextable::bold(ft, part = "header")
ft <- flextable::fontsize(ft, part = "all", size = 10)

if(any(tab_resumo$sig != "")){
  ft <- flextable::add_footer_lines(
    ft, values = c("Significance: * p<0.05; ** p<0.01; *** p<0.001")
  )
  ft <- flextable::align(ft, part = "footer", align = "left")
}

# ---- Exportar para DOCX ----
titulo    <- "Interrupted Time Series (ITS) – Effects of COVID-19 on TB notification rates in Brazil"
subtitulo <- "Estimates from log-scale models (interpreted as percentage changes, 2 decimals)"

doc <- officer::read_docx()
doc <- officer::body_add_par(doc, titulo, style = "heading 1")
doc <- officer::body_add_par(doc, subtitulo, style = "Normal")
doc <- officer::body_add_par(doc, "")
doc <- flextable::body_add_flextable(doc, value = ft)

nota <- paste0(
  "Notes: Estimates on log-scale are expressed as percentage changes: ",
  "Estimate = (exp(coef) − 1) × 100. ‘Immediate level change (2020)’ is the shift at the interruption; ",
  "‘Trend change after 2020’ is the change in post-2020 slope versus pre-2020. ",
  "p-values and percentages are shown with two decimals."
)
doc <- officer::body_add_par(doc, "")
doc <- officer::body_add_par(doc, "Notes", style = "heading 2")
doc <- officer::body_add_par(doc, nota, style = "Normal")

print(doc, target = arq_docx)
message("DOCX salvo em: ", normalizePath(arq_docx))
