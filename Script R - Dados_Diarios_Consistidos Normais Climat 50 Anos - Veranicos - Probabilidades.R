# ============================ PACOTES ============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(broom)
  library(scales)
  library(survival)     # Kaplan–Meier
  # install.packages("evd"); install.packages("ismev"); install.packages("extRemes")
  # library(evd); library(ismev); library(extRemes)   # para GPD (opcional)
  # install.packages("depmixS4")                      # HMM (opcional)
})

# ============================ ENTRADA ============================
arquivo <- "C:/Users/juaci/Downloads/Ocorrencia de Veranicos.txt"

limiar_seco_mm <- 0.5
dur_min_veranico <- 5     # k mínimo (pode variar por análise)
# janelas sazonais úteis (exemplo para safra/estação chuvosa do Cerrado)
meses_janela <- 10:3      # de Out(10) a Mar(3) – ajusta abaixo

# ============================ LEITURA ============================
df <- readr::read_delim(
  arquivo, delim = "\t",
  locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE, trim_ws = TRUE
) |>
  clean_names() |>
  rename(
    mes  = any_of(c("mes","mês")),
    ano  = any_of("ano"),
    data = any_of("data"),
    precipitacao = any_of(c("precipitacao","precipitacao_","precipitação","precipitação_")),
    num_ocorrencia_corrig = any_of("num_ocorrencia_corrig")
  ) |>
  mutate(
    data = suppressWarnings(dmy(data)),
    precipitacao = suppressWarnings(as.numeric(precipitacao))
  ) |>
  filter(!is.na(data)) |>
  arrange(data)

# Ajusta fator mês (1–12)
if(!"mes" %in% names(df)) df <- df |> mutate(mes = month(data))
df <- df |> mutate(ano = if_else(is.na(ano), year(data), as.integer(ano)))

# Marca dia seco
df <- df |>
  mutate(seco = if_else(!is.na(precipitacao) & precipitacao <= limiar_seco_mm, 1L, 0L))

# ====================== DETECÇÃO DE VERANICOS (RUNS) =============
# Função para extrair runs de dias secos (>= dur_min)
runs_veranico <- function(dates, seco_vec, dur_min = 5) {
  stopifnot(length(dates) == length(seco_vec))
  rle_obj <- rle(seco_vec)
  ends <- cumsum(rle_obj$lengths)
  starts <- ends - rle_obj$lengths + 1

  tibble(
    start_id = starts,
    end_id   = ends,
    seco_val = rle_obj$values,
    dur      = rle_obj$lengths
  ) |>
    filter(seco_val == 1, dur >= dur_min) |>
    mutate(
      inicio = dates[start_id],
      fim    = dates[end_id],
      ano_ini = year(inicio),
      mes_ini = month(inicio)
    ) |>
    select(inicio, fim, dur, ano_ini, mes_ini)
}

ver_tbl <- runs_veranico(df$data, df$seco, dur_min = dur_min_veranico)

# =================== PROB. EMPÍRICAS (ANO / MÊS / JANELA) =======
# 1) Prob. anual de ocorrer >=1 veranico com D>=k
prob_anual <- ver_tbl |>
  group_by(ano_ini) |>
  summarise(eventos = n()) |>
  right_join(df |> distinct(ano) |> rename(ano_ini = ano), by = "ano_ini") |>
  mutate(eventos = replace_na(eventos, 0L),
         ocorreu = as.integer(eventos > 0)) |>
  arrange(ano_ini)

p_hat_anual <- mean(prob_anual$ocorreu, na.rm = TRUE)
N_anual <- sum(!is.na(prob_anual$ocorreu))
# IC binomial exato (Clopper–Pearson)
ic_anual <- binom.test(sum(prob_anual$ocorreu, na.rm = TRUE), N_anual)$conf.int
T_ret_anual <- 1 / p_hat_anual

# 2) Prob. mensal: em quantos anos, aquele mês teve >=1 veranico iniciando nele
prob_mensal <- ver_tbl |>
  count(ano_ini, mes_ini, name = "eventos_mes") |>
  mutate(ocorreu = as.integer(eventos_mes > 0)) |>
  group_by(mes_ini) |>
  summarise(
    p_hat = mean(ocorreu, na.rm = TRUE),
    N = n()
  ) |>
  rowwise() |>
  mutate(
    ci = list(binom.test(round(p_hat*N), N)$conf.int),
    p_low = ci[[1]], p_high = ci[[2]],
    retorno_anos = 1/p_hat
  ) |>
  ungroup() |>
  mutate(mes_lab = factor(mes_ini, levels = 1:12, labels = month.abb))

# 3) Prob. sazonal (janela Out–Mar, por exemplo)
is_in_janela <- function(m) {
  # meses_janela = 10:3 (atravessa ano); trata logicamente:
  if (min(meses_janela) <= max(meses_janela)) {
    m %in% meses_janela
  } else {
    m %in% c(meses_janela, meses_janela + 12) # fallback, não usado aqui
  }
}
ver_tbl <- ver_tbl |>
  mutate(na_janela = mes_ini %in% meses_janela)

prob_sazonal <- ver_tbl |>
  filter(na_janela) |>
  count(ano_ini, name = "evt") |>
  mutate(ocorreu = as.integer(evt > 0)) |>
  right_join(df |> distinct(ano) |> rename(ano_ini = ano), by = "ano_ini") |>
  mutate(ocorreu = replace_na(ocorreu, 0L)) |>
  summarise(
    p_hat = mean(ocorreu),
    N = n()
  )
ic_sazonal <- binom.test(round(prob_sazonal$p_hat*prob_sazonal$N), prob_sazonal$N)$conf.int
T_ret_sazonal <- 1 / prob_sazonal$p_hat

# =================== SOBREVIVÊNCIA DAS DURAÇÕES S(d) =============
# Kaplan–Meier para D (todas as durações observadas)
duracoes <- ver_tbl$dur
# Evento sempre "ocorreu" (estamos medindo duração completa), cens = 0
fit_km <- survfit(Surv(time = duracoes, event = rep(1, length(duracoes))) ~ 1, type = "kaplan-meier")

# Função auxiliar: prob de ultrapassar d dias, S(d)
S_empirico <- function(d){
  # KM é por tempo contínuo; durações são discretas, mas aproxima bem
  # Retorna S(d) ~ P(D >= d)
  # aproximação: encontrar S(t) em t = d (ou o mais próximo à direita)
  s <- summary(fit_km, times = d)$surv
  if (length(s) == 0) {
    # se d > max(dur), S = 0
    return(0)
  } else return(s)
}

# Exemplo: prob(D >= 7) e prob(D >= 10)
p_D_ge_7  <- S_empirico(7)
p_D_ge_10 <- S_empirico(10)

# ============ EXTREMOS (POT/GPD) – opcional ======================
# Se quiser modelar cauda das durações (excessos acima u):
# u <- 7
# excessos <- (duracoes[duracoes > u] - u)
# library(evd)
# fit_gpd <- fpot(excessos + u, threshold = u)    # alternativa com "evd"
# # A partir de (xi, beta) obter P(D >= d) na cauda e níveis de retorno:
# # Fórmulas padrão GPD: P(D >= d | D > u) = (1 + xi * (d - u)/beta)^(-1/xi), xi != 0
# # P(D >= d) = P(D > u) * P(D >= d | D > u)

# =================== CONTAGEM ANUAL (>= k) =======================
# nº de eventos por ano
contagem_anual <- ver_tbl |>
  group_by(ano_ini) |>
  summarise(n_eventos = n()) |>
  ungroup()

# Ajustes simples (Poisson/NegBin) – escolha conforme overdispersão:
# library(MASS)
# pois_mod <- glm(n_eventos ~ 1, family = poisson(), data = contagem_anual)
# nb_mod   <- MASS::glm.nb(n_eventos ~ 1, data = contagem_anual)
# # Prob. preditiva de ter >= m eventos:
# m <- 2
# lam <- exp(coef(pois_mod)[1])
# p_ge_m_pois <- 1 - ppois(m-1, lam)

# =================== HMM (wet/dry) – opcional ====================
# # Modela estados latentes (seco/chuvoso) e transições sazonais.
# library(depmixS4)
# hmm_df <- df |>
#   mutate(obs = if_else(precipitacao <= limiar_seco_mm, 1, 0)) |>
#   select(obs)
# set.seed(123)
# mod <- depmix(obs ~ 1, family = binomial(), nstates = 2, data = hmm_df)
# hmm_fit <- fit(mod)
# # A partir de P(seco->seco) = p_DD, a duração de runs secos é ~ geométrica:
# # P(D >= k) = p_DD^(k-1). Para sazonalidade: ajustar por mês (covariáveis) e/ou simular.

# ========================== SAÍDAS-CHAVE =========================
message("=== Probabilidade anual (>= 1 veranico com D >= ", dur_min_veranico, " dias) ===")
print(tibble(
  p_hat = p_hat_anual,
  ic_low = ic_anual[1], ic_high = ic_anual[2],
  retorno_anos = T_ret_anual
))

message("=== Probabilidade mensal (evento iniciando no mês) ===")
prob_mensal_out <- prob_mensal |>
  transmute(
    mes = mes_lab,
    p_hat = p_hat,
    ic_low = p_low, ic_high = p_high,
    retorno_anos
  )
print(prob_mensal_out)

message("=== Probabilidade sazonal (janela) ===")
print(tibble(
  p_hat = prob_sazonal$p_hat,
  ic_low = ic_sazonal[1], ic_high = ic_sazonal[2],
  retorno_anos = T_ret_sazonal
))

message("=== Sobrevivência S(d) = P(D >= d) ===")
print(tibble(
  d = c(5,7,10,15),
  S_d = sapply(c(5,7,10,15), S_empirico)
))





#####################################################################################################





# ====================== PACOTES PARA GRÁFICOS ====================
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(broom)
  library(rlang)
})

# ====================== THEME E DIRETÓRIO ========================
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9, hjust = 0),
    axis.title = element_text(size = 11)
  )

dir.create("figuras", showWarnings = FALSE)

# ======================= 1) CURVA DE SOBREVIVÊNCIA ===============
# KM já ajustado: fit_km (do script anterior)
# Tabela detalhada com IC
km_sum <- summary(fit_km)
km_df <- tibble(
  tempo = km_sum$time,
  S = km_sum$surv,
  S_low = km_sum$lower,
  S_high = km_sum$upper
)

pS <- ggplot(km_df, aes(x = tempo, y = S)) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.2) +
  geom_step(linewidth = 1) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(
    title = "Sobrevivência das durações de veranicos",
    subtitle = "S(d) = Probabilidade de um veranico durar pelo menos d dias (Kaplan–Meier, com IC)",
    x = "Duração (dias)",
    y = "S(d) = P(D ≥ d)",
    caption = "Assume-se que a duração é registrada no dia de encerramento de cada evento."
  ) +
  theme_paper
print(pS)
ggsave("figuras/prob_Sobrevivencia_KM.png", pS, width = 7, height = 4.2, dpi = 300)

# ======================= 2) PROBABILIDADE MENSAL =================
# prob_mensal_out vem do script anterior (com p_hat, IC e retorno_anos)
# Se não existir, recrie a partir de prob_mensal
if (!exists("prob_mensal_out")) {
  prob_mensal_out <- prob_mensal |>
    transmute(
      mes = factor(mes_ini, levels = 1:12, labels = month.abb),
      p_hat = p_hat,
      ic_low = p_low, ic_high = p_high,
      retorno_anos
    )
}

pMes <- ggplot(prob_mensal_out, aes(x = mes, y = p_hat)) +
  geom_col(width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin = ic_low, ymax = ic_high), width = 0.2) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(
    title = "Probabilidade mensal de início de veranico",
    subtitle = "Frequência histórica (com IC binomial exato, Clopper–Pearson)",
    x = "Mês", y = "Probabilidade")

print(pMes)

   
##################################################################################################

# ============================ PACOTES ============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(scales)
  library(broom)
  library(survival)   # Kaplan–Meier
  # ---- opcionais para extremos/HMM ----
  # install.packages(c("evd","ismev","extRemes","depmixS4"))
  # library(evd); library(ismev); library(extRemes); library(depmixS4)
})

# ====================== CONFIGURAÇÕES GERAIS =====================
arquivo <- "C:/Users/juaci/Downloads/Ocorrencia de Veranicos.txt"

limiar_seco_mm     <- 1.0
dur_min_veranico   <- 5      # k mínimo para definir veranico
meses_janela       <- 10:3   # estação chuvosa (Out–Mar). Usado só no rótulo.
use_extremos       <- FALSE  # habilitar GPD (Fig D) se TRUE
use_hmm            <- FALSE  # habilitar HMM (Fig E) se TRUE

# ============================ LEITURA ============================
df <- readr::read_delim(
  arquivo, delim = "\t",
  locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE, trim_ws = TRUE
) |>
  clean_names() |>
  rename(
    mes  = any_of(c("mes","mês")),
    ano  = any_of("ano"),
    data = any_of("data"),
    precipitacao = any_of(c("precipitacao","precipitacao_","precipitação","precipitação_")),
    num_ocorrencia_corrig = any_of("num_ocorrencia_corrig")
  ) |>
  mutate(
    data = suppressWarnings(dmy(data)),
    precipitacao = suppressWarnings(as.numeric(precipitacao))
  ) |>
  filter(!is.na(data)) |>
  arrange(data)

if(!"mes" %in% names(df)) df <- df |> mutate(mes = month(data))
df <- df |> mutate(ano = if_else(is.na(ano), year(data), as.integer(ano)))

# Marca dia seco
df <- df |>
  mutate(seco = if_else(!is.na(precipitacao) & precipitacao <= limiar_seco_mm, 1L, 0L))

# ================== DETECÇÃO DE VERANICOS (RUNS) =================
runs_veranico <- function(dates, seco_vec, dur_min = 5) {
  stopifnot(length(dates) == length(seco_vec))
  rle_obj <- rle(seco_vec)
  ends <- cumsum(rle_obj$lengths)
  starts <- ends - rle_obj$lengths + 1

  tibble(
    start_id = starts,
    end_id   = ends,
    seco_val = rle_obj$values,
    dur      = rle_obj$lengths
  ) |>
    filter(seco_val == 1, dur >= dur_min) |>
    mutate(
      inicio = dates[start_id],
      fim    = dates[end_id],
      ano_ini = year(inicio),
      mes_ini = month(inicio)
    ) |>
    select(inicio, fim, dur, ano_ini, mes_ini)
}
ver_tbl <- runs_veranico(df$data, df$seco, dur_min = dur_min_veranico)
print(ver_tbl)

# ================== PROBABILIDADES (EMPÍRICAS) ===================
# Anual: prob de >= 1 evento no ano
prob_anual <- ver_tbl |>
  count(ano_ini, name = "eventos") |>
  right_join(df |> distinct(ano) |> rename(ano_ini = ano), by = "ano_ini") |>
  mutate(eventos = replace_na(eventos, 0L), ocorreu = as.integer(eventos > 0)) |>
  arrange(ano_ini)
p_hat_anual <- mean(prob_anual$ocorreu, na.rm = TRUE)
N_anual <- sum(!is.na(prob_anual$ocorreu))
ic_anual <- binom.test(sum(prob_anual$ocorreu, na.rm = TRUE), N_anual)$conf.int
T_ret_anual <- 1/p_hat_anual

# Mensal: prob de iniciar evento em cada mês (em anos da série)
prob_mensal <- ver_tbl |>
  count(ano_ini, mes_ini, name = "eventos_mes") |>
  mutate(ocorreu = as.integer(eventos_mes > 0)) |>
  group_by(mes_ini) |>
  summarise(
    p_hat = mean(ocorreu, na.rm = TRUE),
    N = n()
  ) |>
  rowwise() |>
  mutate(
    ci = list(binom.test(round(p_hat*N), N)$conf.int),
    p_low = ci[[1]], p_high = ci[[2]],
    retorno_anos = 1/p_hat
  ) |>
  ungroup() |>
  mutate(mes_lab = factor(mes_ini, levels = 1:12, labels = month.abb))

# Sazonal (Out–Mar)
janela_label <- "Out–Mar"
prob_sazonal <- ver_tbl |>
  mutate(na_janela = mes_ini %in% c(10,11,12,1,2,3)) |>
  filter(na_janela) |>
  count(ano_ini, name = "evt") |>
  mutate(ocorreu = as.integer(evt > 0)) |>
  right_join(df |> distinct(ano) |> rename(ano_ini = ano), by = "ano_ini") |>
  mutate(ocorreu = replace_na(ocorreu, 0L)) |>
  summarise(p_hat = mean(ocorreu), N = n())
ic_sazonal <- binom.test(round(prob_sazonal$p_hat*prob_sazonal$N), prob_sazonal$N)$conf.int
T_ret_sazonal <- 1/prob_sazonal$p_hat

# ==================== SOBREVIVÊNCIA (KAPLAN–MEIER) ===============
duracoes <- ver_tbl$dur
fit_km <- survfit(Surv(time = duracoes, event = rep(1, length(duracoes))) ~ 1, type = "kaplan-meier")
km_df <- broom::tidy(fit_km) |> # time, n.risk, n.event, n.censor, estimate, std.error, conf.low, conf.high
  transmute(
    d = time,
    S = estimate, S_low = conf.low, S_high = conf.high
  )

# ====================== ÍNDICES MENSAIS (HEATMAP) ================
indices_mensais <- df |>
  group_by(ano, mes) |>
  summarise(
    # soma do "tamanho de runs" registrados no mês (aprox intensidade)
    veranicos_dias_total = sum(num_ocorrencia_corrig, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(mes_lab = factor(mes, levels = 1:12, labels = month.abb)) |>
  arrange(ano, mes)

# ========================== ESTILO GRÁFICO =======================
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9),
    axis.title = element_text(size = 11)
  )
dir.create("figuras", showWarnings = FALSE)

# ============================ FIG A ==============================
# Probabilidade mensal de iniciar veranico (≥ k) com IC binomial
figA <- ggplot(prob_mensal, aes(x = mes_lab, y = p_hat)) +
  geom_col(width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin = p_low, ymax = p_high), width = 0.15) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(
    title = paste0("Probabilidade mensal de iniciar veranico (≥ ", dur_min_veranico, " dias)"),
    subtitle = "Estimativa empírica por mês (IC binomial de Clopper–Pearson)",
    x = "Mês", y = "Probabilidade",
    caption = "Probabilidade de observar ≥1 início de veranico no mês (em anos da série)."
  ) +
  theme_paper

figA 

ggsave("figuras/FigA_prob_mensal_IC.png", figA, width = 7, height = 4.2, dpi = 300)

# ============================ FIG B ==============================
# Sobrevivência S(d) = P(D ≥ d) (Kaplan–Meier) com IC
figB <- ggplot(km_df, aes(x = d, y = S)) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.25) +
  geom_step(linewidth = 1) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Sobrevivência das durações de veranico",
    subtitle = "Kaplan–Meier para P(D ≥ d) com faixa de confiança (95%)",
    x = "Duração (dias)", y = "S(d) = Probabilidade acumulada de sobrevivência",
    caption = "Eventos completos (sem censura) construídos pelos runs de dias secos."
  ) +
  theme_paper

figB

ggsave("figuras/FigB_sobrevivencia_KM.png", figB, width = 7, height = 4.2, dpi = 300)

# ============================ FIG C ==============================
# Heatmap ano × mês para “dias em veranicos” (intensidade)
figC <- ggplot(indices_mensais, aes(x = mes_lab, y = ano, fill = veranicos_dias_total)) +
  geom_tile() +
  scale_fill_continuous(type = "viridis", name = "Dias em veranicos") +
  labs(
    title = "Intensidade mensal de veranicos por ano",
    subtitle = "Soma mensal do indicador de duração de veranicos",
    x = "Mês", y = "Ano"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 11)
  )

figC

ggsave("figuras/FigC_heatmap_ano_mes.png", figC, width = 7, height = 6, dpi = 300)

# ============================ FIG D (OPCIONAL) ===================
if (use_extremos) {
  # Exemplo com "extRemes" (ajuste conforme pacote escolhido)
  # library(extRemes)
  # u <- 7
  # exc <- duracoes[duracoes > u] - u
  # if (length(exc) >= 20) {  # mínimo razoável
  #   fit <- extRemes::fevd(exc, type = "GP")
  #   # Curva de nível/retorno (ex.: período em anos vs. duração associada)
  #   # Aqui geramos um gráfico simples apenas ilustrativo:
  #   rp <- c(2,5,10,20,50) # períodos de retorno (anos)
  #   # Converter RP -> quantil de duração (depende da taxa de excedência anual)
  #   # Requer estimar taxa λ de eventos > u por ano. Aproximação:
  #   taxa_anual <- mean(prob_anual$eventos > 0) # eventos/ano (≥ um)
  #   # Este passo pode ser refinado (contagem média de excedências por ano).
  #   # Para ilustração, mantemos a estrutura do gráfico.
  #   figD <- ggplot(data.frame(RP = rp, Duracao = NA), aes(RP, Duracao)) +
  #     geom_line() + geom_point() +
  #     scale_x_log10() +
  #     labs(
  #       title = "Nível de retorno (cauda GPD)",
  #       subtitle = "Exemplo ilustrativo — ajustar conforme a taxa de excedência",
  #       x = "Período de retorno (anos, escala log)", y = "Duração (dias)"
  #     ) + theme_paper
  #   ggsave("figuras/FigD_extremos_GPD.png", figD, width = 7, height = 4.2, dpi = 300)
  # }
}

# ============================ FIG E (OPCIONAL) ===================
if (use_hmm) {
  # # Probabilidade condicional via HMM wet-dry (ilustrativo)
  # library(depmixS4)
  # hmm_df <- df |> transmute(obs = if_else(precipitacao <= limiar_seco_mm, 1, 0))
  # set.seed(123)
  # mod <- depmix(obs ~ 1, family = binomial(), nstates = 2, data = hmm_df)
  # hmm_fit <- fit(mod)
  # # Extrair p_DD (prob de permanecer seco) do estado “seco”
  # # (em modelos completos, usar decoding + sazonalidade por mês)
  # tr <- matrix(getpars(hmm_fit)[grep("trans", names(getpars(hmm_fit)))], nrow = 2, byrow = TRUE)
  # p_DD <- tr[2,2]  # supondo estado 2 = seco (verificar!)
  # # P(D ≥ k | já estamos secos há r) ≈ p_DD^(k - r)
  # r <- 3; k <- 7
  # p_cond <- p_DD^(k - r)
  # figE <- ggplot(data.frame(k = 1:20, p = p_DD^(1:20 - r)), aes(k, p)) +
  #   geom_line() + geom_point() +
  #   scale_y_continuous(labels = label_percent()) +
  #   labs(
  #     title = "Probabilidade condicional de atingir D ≥ k (HMM)",
  #     subtitle = paste0("Condicionado a já estarmos secos há r = ", r, " dias"),
  #     x = "k (dias)", y = "Probabilidade condicional"
  #   ) + theme_paper
  # ggsave("figuras/FigE_prob_condicional_HMM.png", figE, width = 7, height = 4.2, dpi = 300)
}

# ========================== SAÍDAS-RESUMO ========================
cat("\n=== Probabilidade anual (≥ 1 veranico com D ≥", dur_min_veranico, "dias) ===\n")
print(tibble(
  p_hat = p_hat_anual,
  ic_low = ic_anual[1], ic_high = ic_anual[2],
  retorno_anos = T_ret_anual
))

cat("\n=== Probabilidade sazonal ", janela_label, " ===\n", sep = "")
print(tibble(
  p_hat = prob_sazonal$p_hat,
  ic_low = ic_sazonal[1], ic_high = ic_sazonal[2],
  retorno_anos = T_ret_sazonal
))

cat("\nArquivos salvos em: figuras/\n",
    " - FigA_prob_mensal_IC.png\n",
    " - FigB_sobrevivencia_KM.png\n",
    " - FigC_heatmap_ano_mes.png\n",
    if (use_extremos) " - FigD_extremos_GPD.png\n" else "",
    if (use_hmm)     " - FigE_prob_condicional_HMM.png\n" else "",
    sep = "")




##################################################################################################################


# ============================ PACOTES ============================

install.packages("Kendall")
install.packages("trend")

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(scales)
  library(Kendall)    # Mann–Kendall
  library(trend)      # Sen's slope (sens.slope)
  # install.packages("zyp") # alternativa para Sen (zyp.sen), se desejar
})

# ============================ ENTRADA ============================
arquivo <- "C:/Users/juaci/Downloads/Ocorrencia de Veranicos.txt"

limiar_seco_mm   <- 1.0
dur_min_veranico <- 5
salvar_figuras   <- TRUE

# ============================ LEITURA ============================
df <- readr::read_delim(
  arquivo, delim = "\t",
  locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE, trim_ws = TRUE
) |>
  clean_names() |>
  rename(
    mes  = any_of(c("mes","mês")),
    ano  = any_of("ano"),
    data = any_of("data"),
    precipitacao = any_of(c("precipitacao","precipitacao_","precipitação","precipitação_")),
    num_ocorrencia_corrig = any_of("num_ocorrencia_corrig")
  ) |>
  mutate(
    data = suppressWarnings(dmy(data)),
    precipitacao = suppressWarnings(as.numeric(precipitacao))
  ) |>
  filter(!is.na(data)) |>
  arrange(data)

if(!"mes" %in% names(df)) df <- df |> mutate(mes = month(data))
df <- df |> mutate(ano = if_else(is.na(ano), year(data), as.integer(ano)))

# Dia seco (informativo)
df <- df |>
  mutate(seco = if_else(!is.na(precipitacao) & precipitacao <= limiar_seco_mm, 1L, 0L))

# ====================== ÍNDICES ANUAIS/MENSAIS ===================
indices_anuais <- df |>
  group_by(ano) |>
  summarise(
    veranicos_dias_total = sum(num_ocorrencia_corrig, na.rm = TRUE),
    veranico_duracao_max = suppressWarnings(max(num_ocorrencia_corrig, na.rm = TRUE)),
    veranicos_eventos    = sum(num_ocorrencia_corrig > 0, na.rm = TRUE)
  ) |>
  ungroup() |>
  arrange(ano)

indices_mensais <- df |>
  group_by(ano, mes) |>
  summarise(
    veranicos_dias_total = sum(num_ocorrencia_corrig, na.rm = TRUE),
    veranico_duracao_max = suppressWarnings(max(num_ocorrencia_corrig, na.rm = TRUE)),
    veranicos_eventos    = sum(num_ocorrencia_corrig > 0, na.rm = TRUE)
  ) |>
  ungroup() |>
  arrange(ano, mes)

# ===================== FUNÇÕES DE TENDÊNCIA ======================
mk_sen <- function(x, y) {
  # Remove NA e garante ordenação por x
  df2 <- tibble(x = x, y = y) |> drop_na() |> arrange(x)
  if (nrow(df2) < 8) return(tibble(
    tau = NA_real_, p_value = NA_real_,
    sen_slope = NA_real_, sen_ci_low = NA_real_, sen_ci_high = NA_real_
  ))
  mk <- Kendall::MannKendall(df2$y)
  ss <- tryCatch(trend::sens.slope(df2$y ~ df2$x), error = function(e) NULL)
  if (is.null(ss)) {
    tibble(tau = mk$tau, p_value = mk$sl,
           sen_slope = NA_real_, sen_ci_low = NA_real_, sen_ci_high = NA_real_)
  } else {
    tibble(tau = mk$tau, p_value = mk$sl,
           sen_slope = unname(ss$estimates[1]),
           sen_ci_low = unname(ss$conf.int[1]),
           sen_ci_high = unname(ss$conf.int[2]))
  }
}

mk_sen_by_month <- function(df_m, var) {
  df_m |>
    group_by(mes) |>
    group_modify(~{
      res <- mk_sen(.x$ano, .x[[var]])
      res
    }) |>
    ungroup() |>
    mutate(var = var)
}

# ====================== TESTES ANUAIS (MK + SEN) =================
vars <- c("veranicos_dias_total", "veranico_duracao_max", "veranicos_eventos")

trend_anuais <- map_dfr(vars, ~{
  res <- mk_sen(indices_anuais$ano, indices_anuais[[.x]])
  res |> mutate(var = .x)
}) |>
  relocate(var)

print(trend_anuais)
# Interpretação: tau>0 tendência crescente; p_value<0.05 significativo.
# sen_slope ~ variação média por ano (escala da variável).

# ==================== TESTES MENSAIS (MK + SEN) ==================
trend_mensais <- map_dfr(vars, ~ mk_sen_by_month(indices_mensais, .x)) |>
  relocate(var, mes)

print(trend_mensais)
# Interpretação por mês (ao longo dos anos para cada mês específico).

# =================== LOESS + TENDÊNCIA POR DÉCADA ================
# LOESS para visualização + regressão linear por década (inclinação/ano)
indices_anuais <- indices_anuais |>
  mutate(decada = paste0(floor(ano/10)*10, "s"))

# Função para regressão por década
lm_by_decade <- function(dat, var) {
  dat |>
    group_by(decada) |>
    summarise(
      n_anos = n(),
      slope_ano = ifelse(n_anos >= 5,
                         coef(lm(dat[[var]][dat$decada == unique(decada)] ~
                                 dat$ano[dat$decada == unique(decada)]))[2],
                         NA_real_),
      p_value = ifelse(n_anos >= 5,
                       summary(lm(dat[[var]][dat$decada == unique(decada)] ~
                                  dat$ano[dat$decada == unique(decada)]))$coefficients[2,4],
                       NA_real_)
    ) |>
    ungroup() |>
    mutate(var = var, .before = 1)
}

decadal_slopes <- map_dfr(vars, ~ lm_by_decade(indices_anuais, .x))
print(decadal_slopes)

# ============================ GRÁFICOS ============================
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9),
    axis.title = element_text(size = 11)
  )

plot_loess <- function(dat, yvar, ylab) {
  ggplot(dat, aes(x = ano, y = .data[[yvar]])) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1.4, alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, span = 0.3, linewidth = 0.9) +
    labs(
      title = paste0("Série anual — ", ylab),
      subtitle = "Suavização LOESS (faixa = IC aproximado)",
      x = "Ano", y = ylab,
      caption = "Teste de tendência reportado por Mann–Kendall e Sen's slope no texto/tabela."
    ) +
    theme_paper
}

p1 <- plot_loess(indices_anuais, "veranicos_dias_total",
                 "Dias em veranicos (total no ano)")
p2 <- plot_loess(indices_anuais, "veranico_duracao_max",
                 "Duração máxima de veranico (dias)")
p3 <- plot_loess(indices_anuais, "veranicos_eventos",
                 "Número de eventos de veranico (n/ano)")

if (salvar_figuras) {
  dir.create("figuras", showWarnings = FALSE)
  ggsave("figuras/TrendA_loess_dias_total.png", p1, width = 7, height = 4.2, dpi = 300)
  ggsave("figuras/TrendB_loess_duracao_max.png", p2, width = 7, height = 4.2, dpi = 300)
  ggsave("figuras/TrendC_loess_eventos.png", p3, width = 7, height = 4.2, dpi = 300)
}

# ==================== SAÍDA RESUMIDA PARA TEXTO ==================
cat("\n==== ANUAL: Mann–Kendall + Sen's slope ====\n")
trend_anuais |>
  mutate(
    direcao = case_when(
      is.na(tau) ~ NA_character_,
      tau > 0 ~ "crescente",
      tau < 0 ~ "decrescente",
      TRUE ~ "nula"
    ),
    signif = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.05 ~ "significativa (p<0.05)",
      p_value < 0.10 ~ "marginal (0.05≤p<0.10)",
      TRUE ~ "não significativa"
    )
  ) |>
  print(n = Inf)

cat("\n==== MENSAL: Mann–Kendall + Sen's slope (por mês) ====\n")
trend_mensais |>
  mutate(
    mes_lab = factor(mes, levels = 1:12, labels = month.abb),
    direcao = case_when(
      is.na(tau) ~ NA_character_,
      tau > 0 ~ "crescente",
      tau < 0 ~ "decrescente",
      TRUE ~ "nula"
    ),
    signif = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.05 ~ "significativa (p<0.05)",
      p_value < 0.10 ~ "marginal (0.05≤p<0.10)",
      TRUE ~ "não significativa"
    )
  ) |>
  arrange(var, mes) |>
  print(n = Inf)

cat("\n==== DÉCADAS: inclinação/ano por década (LM) ====\n")
decadal_slopes |>
  mutate(
    signif = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.05 ~ "significativa (p<0.05)",
      p_value < 0.10 ~ "marginal (0.05≤p<0.10)",
      TRUE ~ "não significativa"
    )
  ) |>
  arrange(var, decada) |>
  print(n = Inf)



###############################################################################################################

# ============================ PACOTES ============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(scales)
})

# ========================== PARÂMETROS ===========================
arquivo <- "C:/Users/juaci/Downloads/Ocorrencia de Veranicos.txt"


limiar_seco_mm     <- 1.0
dur_min_veranico   <- 5         # duração mínima para ser veranico
mes_alvo           <- 7         # 1=Jan, 2=Fev, ..., 12=Dez
criterio_mes       <- "interseccao"  # "inicio", "fim" ou "interseccao"
salvar_figura      <- TRUE

# ====================== LEITURA E PREPARO ========================
df <- readr::read_delim(
  arquivo, delim = "\t",
  locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE, trim_ws = TRUE
) |>
  clean_names() |>
  rename(
    mes  = any_of(c("mes","mês")),
    ano  = any_of("ano"),
    data = any_of("data"),
    precipitacao = any_of(c("precipitacao","precipitacao_","precipitação","precipitação_")),
    num_ocorrencia_corrig = any_of("num_ocorrencia_corrig")
  ) |>
  mutate(
    data = suppressWarnings(dmy(data)),
    precipitacao = suppressWarnings(as.numeric(precipitacao))
  ) |>
  filter(!is.na(data)) |>
  arrange(data)

if (!"mes" %in% names(df)) df <- df |> mutate(mes = month(data))
df <- df |> mutate(ano = if_else(is.na(ano), year(data), as.integer(ano)))

# Dia seco (definição operacional)
df <- df |>
  mutate(seco = if_else(!is.na(precipitacao) & precipitacao <= limiar_seco_mm, 1L, 0L))

# ============ DETECÇÃO DE VERANICOS (runs de dias secos) =========
runs_veranico <- function(dates, seco_vec, dur_min = 5) {
  stopifnot(length(dates) == length(seco_vec))
  rle_obj <- rle(seco_vec)
  ends   <- cumsum(rle_obj$lengths)
  starts <- ends - rle_obj$lengths + 1

  tibble(
    start_id = starts,
    end_id   = ends,
    seco_val = rle_obj$values,
    dur      = rle_obj$lengths
  ) |>
    filter(seco_val == 1, dur >= dur_min) |>
    mutate(
      inicio  = dates[start_id],
      fim     = dates[end_id],
      ano_ini = year(inicio),
      mes_ini = month(inicio),
      ano_fim = year(fim),
      mes_fim = month(fim)
    ) |>
    select(inicio, fim, dur, ano_ini, mes_ini, ano_fim, mes_fim)
}

ver_tbl <- runs_veranico(df$data, df$seco, dur_min = dur_min_veranico)

# ========== FUNÇÃO: MÁXIMO ANUAL RESTRITO A UM MÊS ===============
max_veranico_por_ano_mes <- function(events_tbl, mes_alvo, criterio = c("inicio","fim","interseccao")) {
  criterio <- match.arg(criterio)

  sel <- switch(
    criterio,
    "inicio" = events_tbl |> filter(mes_ini == mes_alvo) |> mutate(ano_ref = ano_ini),
    "fim"    = events_tbl |> filter(mes_fim == mes_alvo) |> mutate(ano_ref = ano_fim),
    "interseccao" = {
      # considera qualquer run que tenha pelo menos 1 dia dentro do mês-alvo
      # aproximação robusta via sequência de meses entre início e fim
      events_tbl |>
        rowwise() |>
        mutate(intercepta =
                 any(month(seq(inicio, fim, by = "1 day")) == mes_alvo)) |>
        ungroup() |>
        filter(intercepta) |>
        # para o ano de referência, usa o ano do mês-alvo dentro do intervalo
        rowwise() |>
        mutate(
          anos_meses = list(tibble(dt = seq(inicio, fim, by = "1 day"),
                                   a = year(dt), m = month(dt))),
          ano_ref = anos_meses$a[which(anos_meses$m == mes_alvo)][1] # primeiro ano em que toca o mês
        ) |>
        ungroup() |>
        select(-anos_meses)
    }
  )

  sel |>
    group_by(ano_ref) |>
    summarise(veranico_duracao_max_mes = max(dur, na.rm = TRUE), .groups = "drop") |>
    arrange(ano_ref) |>
    rename(ano = ano_ref)
}

max_por_ano_mes <- max_veranico_por_ano_mes(ver_tbl, mes_alvo = mes_alvo, criterio = criterio_mes)

# ========================== VISUALIZAÇÃO ==========================
lab_mes <- month.abb[mes_alvo]  # rótulo (Jan, Feb, ...)

p <- ggplot(max_por_ano_mes, aes(x = ano, y = veranico_duracao_max_mes)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_point(size = 1.6, color = "black") +
  geom_smooth(method = "loess", se = TRUE, span = 0.35, linewidth = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = paste0("Duração máxima de veranicos por ano — mês: ", lab_mes),
    subtitle = paste0("Critério: ", criterio_mes,
                      " | Definição de veranico: runs secos ≥ ", dur_min_veranico, " dias (≤ ",
                      limiar_seco_mm, " mm/dia)"),
    x = "Ano", y = "Duração máxima (dias)",
    caption = "Base diária; runs de dias secos detectados por RLE. Linha azul: LOESS (IC sombreado)."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p)
if (salvar_figura) {
  dir.create("figuras", showWarnings = FALSE)
  ggsave(
    filename = sprintf("figuras/duracao_max_por_ano_mes_%s_%s.png", lab_mes, criterio_mes),
    plot = p, width = 9, height = 3.8, dpi = 300
  )
}

# =================== SAÍDA TABULAR (opcional) ====================
max_por_ano_mes |> head()


#############################################################################################################


# ============================ PACOTES ============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(scales)
  library(Kendall)   # Mann–Kendall
  library(trend)     # Sen's slope (sens.slope)
})

# ======================== PARÂMETROS GERAIS ======================
arquivo <- "C:/Users/juaci/Downloads/Ocorrencia de Veranicos.txt"

limiar_seco_mm    <- 1.0
dur_min_veranico  <- 5          # definição operacional (>=5 dias)
mes_alvo          <- 7          # 1=Jan ... 12=Dez (ex.: 7=Julho)
criterio_mes      <- "interseccao"  # "inicio", "fim" ou "interseccao"
salvar_figuras    <- TRUE

# ============================== LEITURA ==========================
df <- readr::read_delim(
  arquivo, delim = "\t",
  locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE, trim_ws = TRUE
) |>
  clean_names() |>
  rename(
    mes  = any_of(c("mes","mês")),
    ano  = any_of("ano"),
    data = any_of("data"),
    precipitacao = any_of(c("precipitacao","precipitacao_","precipitação","precipitação_")),
    num_ocorrencia_corrig = any_of("num_ocorrencia_corrig")
  ) |>
  mutate(
    data = suppressWarnings(dmy(data)),
    precipitacao = suppressWarnings(as.numeric(precipitacao))
  ) |>
  filter(!is.na(data)) |>
  arrange(data)

if (!"mes" %in% names(df)) df <- df |> mutate(mes = month(data))
df <- df |> mutate(ano = if_else(is.na(ano), year(data), as.integer(ano)))

# Marca dia seco
df <- df |>
  mutate(seco = if_else(!is.na(precipitacao) & precipitacao <= limiar_seco_mm, 1L, 0L))

# ============ DETECÇÃO DE VERANICOS (runs de dias secos) =========
runs_veranico <- function(dates, seco_vec, dur_min = 5) {
  stopifnot(length(dates) == length(seco_vec))
  rle_obj <- rle(seco_vec)
  ends   <- cumsum(rle_obj$lengths)
  starts <- ends - rle_obj$lengths + 1

  tibble(
    start_id = starts,
    end_id   = ends,
    seco_val = rle_obj$values,
    dur      = rle_obj$lengths
  ) |>
    filter(seco_val == 1, dur >= dur_min) |>
    mutate(
      inicio  = dates[start_id],
      fim     = dates[end_id],
      ano_ini = year(inicio),
      mes_ini = month(inicio),
      ano_fim = year(fim),
      mes_fim = month(fim)
    ) |>
    select(inicio, fim, dur, ano_ini, mes_ini, ano_fim, mes_fim)
}

ver_tbl <- runs_veranico(df$data, df$seco, dur_min = dur_min_veranico)

# ====== SÉRIE: DURAÇÃO MÁXIMA POR ANO RESTRITA AO MÊS ESCOLHIDO ==
max_veranico_por_ano_mes <- function(events_tbl, mes_alvo, criterio = c("inicio","fim","interseccao")) {
  criterio <- match.arg(criterio)

  sel <- switch(
    criterio,
    "inicio" = events_tbl |> filter(mes_ini == mes_alvo) |> mutate(ano_ref = ano_ini),
    "fim"    = events_tbl |> filter(mes_fim == mes_alvo) |> mutate(ano_ref = ano_fim),
    "interseccao" = {
      events_tbl |>
        rowwise() |>
        mutate(intercepta = any(month(seq(inicio, fim, by = "1 day")) == mes_alvo)) |>
        ungroup() |>
        filter(intercepta) |>
        # ano de referência = ano dentro do evento quando ele toca o mês-alvo
        rowwise() |>
        mutate(
          anos_meses = list(tibble(dt = seq(inicio, fim, by = "1 day"),
                                   a = year(dt), m = month(dt))),
          ano_ref = anos_meses$a[which(anos_meses$m == mes_alvo)][1]
        ) |>
        ungroup() |>
        select(-anos_meses)
    }
  )

  sel |>
    group_by(ano_ref) |>
    summarise(duracao_max = max(dur, na.rm = TRUE), .groups = "drop") |>
    arrange(ano_ref) |>
    rename(ano = ano_ref)
}

serie_max <- max_veranico_por_ano_mes(ver_tbl, mes_alvo = mes_alvo, criterio = criterio_mes)

# ======================= TESTE MANN–KENDALL ======================
mk_df <- { 
  d <- drop_na(serie_max)
  if (nrow(d) < 8) {
    tibble(tau = NA_real_, p_value = NA_real_)
  } else {
    mk <- Kendall::MannKendall(d$duracao_max)
    tibble(tau = as.numeric(mk$tau), p_value = as.numeric(mk$sl))
  }
}


# ========================= SEN'S SLOPE (robusto) =================
# Requer: serie_max com colunas ano (num/int) e duracao_max (num)
sen_df <- {
  d <- serie_max |>
    dplyr::arrange(ano) |>
    dplyr::select(ano, duracao_max) |>
    tidyr::drop_na()

  if (nrow(d) < 8) {
    tibble::tibble(
      sen_slope  = NA_real_,
      sen_ci_low = NA_real_,
      sen_ci_high= NA_real_,
      metodo     = NA_character_
    )
  } else {
    d$ano         <- as.numeric(d$ano)
    d$duracao_max <- as.numeric(d$duracao_max)

    # 1) Tenta zyp::zyp.sen (preferível; aceita fórmula y ~ x)
    if (requireNamespace("zyp", quietly = TRUE)) {
      fit <- zyp::zyp.sen(duracao_max ~ ano, data = d)
      tibble::tibble(
        sen_slope  = unname(stats::coef(fit)[["ano"]]),
        sen_ci_low = NA_real_,   # zyp.sen não retorna IC direto
        sen_ci_high= NA_real_,
        metodo     = "zyp.sen"
      )
    } else {
      # 2) Fallback: trend::sens.slope sobre um índice sequencial (t)
      # evita o bug do método com fórmula
      d <- d |>
        dplyr::mutate(t = dplyr::row_number())
      ss <- trend::sens.slope(d$duracao_max ~ d$t)
      tibble::tibble(
        sen_slope  = unname(ss$estimates[1]),
        sen_ci_low = unname(ss$conf.int[1]),
        sen_ci_high= unname(ss$conf.int[2]),
        metodo     = "trend::sens.slope (t-index)"
      )
    }
  }
}

# =============== REGRESSÃO POR DÉCADA (SLOPE/ANO) ===============
serie_max <- serie_max |>
  mutate(decada = paste0(floor(ano/10)*10, "s"))

decadal_slopes <- serie_max |>
  group_by(decada) |>
  summarise(
    n_anos = n(),
    slope_ano = ifelse(n_anos >= 5,
                       coef(lm(duracao_max ~ ano))[2], NA_real_),
    p_value = ifelse(n_anos >= 5,
                     summary(lm(duracao_max ~ ano))$coefficients[2,4], NA_real_)
  ) |>
  ungroup()

# ============================= GRÁFICO ===========================
lab_mes <- month.abb[mes_alvo]
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9),
    axis.title = element_text(size = 11)
  )

p <- ggplot(serie_max, aes(x = ano, y = duracao_max)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_point(size = 1.6, color = "black") +
  geom_smooth(method = "loess", se = TRUE, span = 0.35, linewidth = 0.9) +
  labs(
    title = paste0("Duração máxima de veranicos por ano — mês: ", lab_mes),
    subtitle = paste0("Critério: ", criterio_mes,
                      " | Definição: runs secos ≥ ", dur_min_veranico, " dias (≤ ",
                      limiar_seco_mm, " mm/dia)"),
    x = "Ano", y = "Duração máxima (dias)",
    caption = "Tendência formal: Mann–Kendall + Sen's slope; LOESS para visualização."
  ) +
  theme_paper

if (salvar_figuras) {
  dir.create("figuras", showWarnings = FALSE)
  ggsave(sprintf("figuras/MaxVeranico_%s_%s.png", lab_mes, criterio_mes),
         p, width = 9, height = 4, dpi = 300)
}
print(p)


# ============================ SAÍDAS =============================
cat("\n=== Tendência global (duração máxima, mês =", lab_mes, ") ===\n")

# Garante objetos existentes
if (!exists("mk_df") || nrow(mk_df) == 0) {
  mk_df <- tibble::tibble(tau = NA_real_, p_value = NA_real_)
}

dplyr::bind_cols(mk_df, sen_df) |>
  dplyr::mutate(
    direcao = dplyr::case_when(
      is.na(tau) ~ NA_character_,
      tau > 0    ~ "crescente",
      tau < 0    ~ "decrescente",
      TRUE       ~ "nula"
    ),
    signif = dplyr::case_when(
      is.na(p_value)     ~ NA_character_,
      p_value < 0.05     ~ "significativa (p<0.05)",
      p_value < 0.10     ~ "marginal (0.05≤p<0.10)",
      TRUE               ~ "não significativa"
    )
  ) |>
  print(n = Inf)



##################################################################################################
##################################################################################################

# ============================ PACOTES ============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(scales)
  library(Kendall)   # Mann–Kendall
  library(trend)     # Sen's slope (sens.slope)
  # opcional (preferível para Theil–Sen): install.packages("zyp")
  # library(zyp)
})

# ======================== PARÂMETROS GERAIS ======================
arquivo <- "C:/Users/juaci/Downloads/Ocorrencia de Veranicos.txt"

limiar_seco_mm    <- 1.0
dur_min_veranico  <- 5                # definição operacional (>= 5 dias)
mes_alvo          <- 7                # 1=Jan ... 12=Dez (ex.: 7=Julho)
criterio_mes      <- "interseccao"    # "inicio", "fim" ou "interseccao"
salvar_figuras    <- TRUE
B_boot            <- 1000             # bootstrap p/ plano C (Theil–Sen manual)

# ============================== LEITURA ==========================
df <- readr::read_delim(
  arquivo, delim = "\t",
  locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE, trim_ws = TRUE
) |>
  clean_names() |>
  rename(
    mes  = any_of(c("mes","mês")),
    ano  = any_of("ano"),
    data = any_of("data"),
    precipitacao = any_of(c("precipitacao","precipitacao_","precipitação","precipitação_")),
    num_ocorrencia_corrig = any_of("num_ocorrencia_corrig")
  ) |>
  mutate(
    data = suppressWarnings(dmy(data)),
    precipitacao = suppressWarnings(as.numeric(precipitacao))
  ) |>
  filter(!is.na(data)) |>
  arrange(data)

if (!"mes" %in% names(df)) df <- df |> mutate(mes = month(data))
df <- df |> mutate(ano = if_else(is.na(ano), year(data), as.integer(ano)))

# Marca dia seco
df <- df |>
  mutate(seco = if_else(!is.na(precipitacao) & precipitacao <= limiar_seco_mm, 1L, 0L))

# ============ DETECÇÃO DE VERANICOS (runs de dias secos) =========
runs_veranico <- function(dates, seco_vec, dur_min = 5) {
  stopifnot(length(dates) == length(seco_vec))
  rle_obj <- rle(seco_vec)
  ends   <- cumsum(rle_obj$lengths)
  starts <- ends - rle_obj$lengths + 1

  tibble(
    start_id = starts,
    end_id   = ends,
    seco_val = rle_obj$values,
    dur      = rle_obj$lengths
  ) |>
    filter(seco_val == 1, dur >= dur_min) |>
    mutate(
      inicio  = dates[start_id],
      fim     = dates[end_id],
      ano_ini = year(inicio),
      mes_ini = month(inicio),
      ano_fim = year(fim),
      mes_fim = month(fim)
    ) |>
    select(inicio, fim, dur, ano_ini, mes_ini, ano_fim, mes_fim)
}

ver_tbl <- runs_veranico(df$data, df$seco, dur_min = dur_min_veranico)

# ====== SÉRIE: DURAÇÃO MÁXIMA POR ANO RESTRITA AO MÊS ESCOLHIDO ==
max_veranico_por_ano_mes <- function(events_tbl, mes_alvo, criterio = c("inicio","fim","interseccao")) {
  criterio <- match.arg(criterio)

  sel <- switch(
    criterio,
    "inicio" = events_tbl |> filter(mes_ini == mes_alvo) |> mutate(ano_ref = ano_ini),
    "fim"    = events_tbl |> filter(mes_fim == mes_alvo) |> mutate(ano_ref = ano_fim),
    "interseccao" = {
      events_tbl |>
        rowwise() |>
        mutate(intercepta = any(month(seq(inicio, fim, by = "1 day")) == mes_alvo)) |>
        ungroup() |>
        filter(intercepta) |>
        rowwise() |>
        mutate(
          anos_meses = list(tibble(dt = seq(inicio, fim, by = "1 day"),
                                   a = year(dt), m = month(dt))),
          ano_ref = anos_meses$a[which(anos_meses$m == mes_alvo)][1]
        ) |>
        ungroup() |>
        select(-anos_meses)
    }
  )

  sel |>
    group_by(ano_ref) |>
    summarise(duracao_max = max(dur, na.rm = TRUE), .groups = "drop") |>
    arrange(ano_ref) |>
    rename(ano = ano_ref)
}

serie_max <- max_veranico_por_ano_mes(ver_tbl, mes_alvo = mes_alvo, criterio = criterio_mes)

# ======================= TESTE MANN–KENDALL ======================
mk_df <- {
  d <- drop_na(serie_max)
  if (nrow(d) < 8) {
    tibble(tau = NA_real_, p_value = NA_real_)
  } else {
    mk <- Kendall::MannKendall(d$duracao_max)
    tibble(tau = as.numeric(mk$tau), p_value = as.numeric(mk$sl))
  }
}

# ======================= THEIL–SEN (ROBUSTO) =====================
# Plano C: função manual com IC por bootstrap (independe de pacotes)
theilsen_manual <- function(x, y, B = 1000, seed = 123){
  stopifnot(length(x) == length(y))
  o <- order(x); x <- as.numeric(x[o]); y <- as.numeric(y[o])
  n <- length(x)
  if (n < 2) return(list(slope = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  # todos os pares i<j
  S <- c()
  for (i in 1:(n-1)){
    dx <- x[(i+1):n] - x[i]
    dy <- y[(i+1):n] - y[i]
    S <- c(S, dy/dx)
  }
  slope_hat <- median(S, na.rm = TRUE)

  # IC bootstrap (percentil)
  set.seed(seed)
  boot <- numeric(B)
  for (b in 1:B){
    idx <- sample.int(n, n, replace = TRUE)
    xb <- x[idx]; yb <- y[idx]
    Sb <- c()
    for (i in 1:(n-1)){
      dx <- xb[(i+1):n] - xb[i]
      dy <- yb[(i+1):n] - yb[i]
      Sb <- c(Sb, dy/dx)
    }
    boot[b] <- median(Sb, na.rm = TRUE)
  }
  ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  list(slope = slope_hat, ci_low = ci[1], ci_high = ci[2])
}

# Sen's slope com fallbacks: zyp -> trend(data=, t-index) -> manual
sen_df <- {
  d <- serie_max |>
    arrange(ano) |>
    select(ano, duracao_max) |>
    drop_na()

  if (nrow(d) < 8) {
    tibble(sen_slope = NA_real_, sen_ci_low = NA_real_, sen_ci_high = NA_real_, metodo = NA_character_)
  } else {
    d$ano         <- as.numeric(d$ano)
    d$duracao_max <- as.numeric(d$duracao_max)

    if (requireNamespace("zyp", quietly = TRUE)) {
      fit <- zyp::zyp.sen(duracao_max ~ ano, data = d)
      tibble(
        sen_slope  = unname(stats::coef(fit)[["ano"]]),
        sen_ci_low = NA_real_,  # zyp.sen não retorna IC direto
        sen_ci_high= NA_real_,
        metodo     = "zyp.sen"
      )
    } else {
      # tenta trend::sens.slope usando data = d e índice t
      ok <- TRUE
      out <- try({
        d2 <- d |> mutate(t = row_number())
        trend::sens.slope(duracao_max ~ t, data = d2)
      }, silent = TRUE)

      if (inherits(out, "try-error")) ok <- FALSE

      if (ok) {
        tibble(
          sen_slope  = unname(out$estimates[1]),
          sen_ci_low = unname(out$conf.int[1]),
          sen_ci_high= unname(out$conf.int[2]),
          metodo     = "trend::sens.slope (t-index)"
        )
      } else {
        th <- theilsen_manual(d$ano, d$duracao_max, B = B_boot)
        tibble(
          sen_slope  = unname(th$slope),
          sen_ci_low = unname(th$ci_low),
          sen_ci_high= unname(th$ci_high),
          metodo     = "Theil–Sen manual (bootstrap CI)"
        )
      }
    }
  }
}

# =============== REGRESSÃO POR DÉCADA (SLOPE/ANO) ===============
serie_max <- serie_max |>
  mutate(decada = paste0(floor(ano/10)*10, "s"))

decadal_slopes <- serie_max |>
  group_by(decada) |>
  summarise(
    n_anos = n(),
    slope_ano = ifelse(n_anos >= 5, coef(lm(duracao_max ~ ano))[2], NA_real_),
    p_value   = ifelse(n_anos >= 5, summary(lm(duracao_max ~ ano))$coefficients[2,4], NA_real_)
  ) |>
  ungroup()

# ============================= GRÁFICO ===========================
lab_mes <- month.abb[mes_alvo]
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9),
    axis.title = element_text(size = 11)
  )

p <- ggplot(serie_max, aes(x = ano, y = duracao_max)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_point(size = 1.6, color = "black") +
  geom_smooth(method = "loess", se = TRUE, span = 0.35, linewidth = 0.9) +
  labs(
    title = paste0("Duração máxima de veranicos por ano — mês: ", lab_mes),
    subtitle = paste0("Critério: ", criterio_mes,
                      " | Definição: runs secos ≥ ", dur_min_veranico, " dias (≤ ",
                      limiar_seco_mm, " mm/dia)"),
    x = "Ano", y = "Duração máxima (dias)",
    caption = "Tendência formal: Mann–Kendall + Sen's slope; LOESS para visualização."
  ) +
  theme_paper

if (salvar_figuras) {
  dir.create("figuras", showWarnings = FALSE)
  ggsave(sprintf("figuras/MaxVeranico_%s_%s.png", lab_mes, criterio_mes),
         p, width = 9, height = 4, dpi = 300)
}
print(p)

# ============================ SAÍDAS =============================
cat("\n=== Tendência global (duração máxima, mês =", lab_mes, ") ===\n")
if (!exists("mk_df") || nrow(mk_df) == 0) mk_df <- tibble(tau = NA_real_, p_value = NA_real_)

bind_cols(mk_df, sen_df) |>
  mutate(
    direcao = case_when(
      is.na(tau) ~ NA_character_,
      tau > 0    ~ "crescente",
      tau < 0    ~ "decrescente",
      TRUE       ~ "nula"
    ),
    signif = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.05 ~ "significativa (p<0.05)",
      p_value < 0.10 ~ "marginal (0.05≤p<0.10)",
      TRUE           ~ "não significativa"
    )
  ) |>
  print(n = Inf)

cat("\n=== Slope por década (LM) ===\n")
decadal_slopes |>
  mutate(
    signif = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.05 ~ "significativa (p<0.05)",
      p_value < 0.10 ~ "marginal (0.05≤p<0.10)",
      TRUE           ~ "não significativa"
    )
  ) |>
  arrange(decada) |>
  print(n = Inf)


###############################################################################







