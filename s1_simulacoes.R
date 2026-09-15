# Pacotes -----------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  kableExtra,
  foreign,
  tidyverse,
  DescTools,
  questionr,
  stats,
  survival,
  survminer,
  broom.helpers,
  gtsummary,
  gt,
  purrr,
  naniar,
  patchwork,
  mice,
  parallel,
  pbapply
)

# Simulação Missing ------------------------------------------------------
rodar_sim <- function(porc, metodo, dados_base, padrao_mat) {
  simulacao <- ampute(data = dados_base, prop = porc, mech = metodo, patterns = padrao_mat)
  return(simulacao$amp)
}

# Deleção Listwise ----------------------------------------------
ajustar_listwise_unico <- function(df, betas_reais_ref) {
  df_omit <- na.omit(df)
  mod <- tryCatch({
    survival::survreg(
      survival::Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + CATEATEND +
        RECIDIVA + METASTASE + QUIMIO + ESCOLARI + LEI_DIAGTRAT,
      data = df_omit, dist = "loglogis"
    )
  }, error = function(e) return(NA))
  
  beta_res <- if (length(mod) == 1 && is.na(mod)) rep(NA, length(betas_reais_ref)) else coef(mod)
  return(list(df_omit = df_omit,
              modelo = mod, 
              beta = beta_res))
}

# Imputação Múltipla --------------------------------------------
ajustar_mice_unico <- function(df, num_imp) {
  mice_res <- mice(df, m = num_imp, method = "pmm", printFlag = FALSE)
  modelos_brutos <- with(mice_res, 
                         survival::survreg(
                           survival::Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + CATEATEND +
                             RECIDIVA + METASTASE + QUIMIO + ESCOLARI + LEI_DIAGTRAT,
                           dist = "loglogis"
                         ))
  modelo_unificado <- pool(modelos_brutos)
  return(list(dados_imputados = mice_res, 
              modelos_brutos = modelos_brutos,
              modelo_pool = modelo_unificado))
}
# Dados --------------------------------------------------------
dados <- read.csv("dados_prostata_completos.csv")

colunas_categoricas <- sapply(dados, is.character)
dados[colunas_categoricas] <- lapply(dados[colunas_categoricas], as.factor)

# Banco temporário para os dados faltantes
banco_temp <- dados

padrao <- rep(1, ncol(banco_temp))
padrao[colnames(banco_temp) == "IDADE"] <- 0
matriz_padrao <- matrix(padrao, nrow = 1)

# Modelo Completo ---------------------------------------------
modelo_completo <- survreg(
  Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + CATEATEND +
    RECIDIVA + METASTASE + QUIMIO + ESCOLARI + LEI_DIAGTRAT,
  data = dados, 
  dist = "loglogis"
)

betas_reais <- coef(modelo_completo)
betas_vcov <- modelo_completo$var
nomes_variaveis <- names(betas_reais)

df_referencia <- data.frame(
  Variavel = names(betas_reais),
  Beta_Real = betas_reais
)

se_idade_real <- summary(modelo_completo)$table["IDADE", "Std. Error"]

tbl_regression(
  modelo_completo, 
  exponentiate = TRUE,
  estimate_fun = ~ style_number(.x, digits = 4),
) |> 
  bold_p(t = 0.05) |> 
  bold_labels()

# Cenários
cenarios <- expand.grid(
  porc = c(0.05, 0.1, 0.2, 0.4, 0.6), # ALTERAR PARA CADA COMPUTADOR
  mech = c("MAR"), # ALTERAR PARA CADA COMPUTADOR
  stringsAsFactors = FALSE
)

# Tabela para guardar os tempos computacionais
registro_tempos <- cenarios
registro_tempos$tempo_amputacao <- NA
registro_tempos$tempo_listwise <- NA
registro_tempos$tempo_mice <- NA

if(!dir.exists("Resultados")) dir.create("Resultados")

# Simulando dados --------------------------------------------
# Amputação
cat("Iniciando Fase 1: Amputação...\n")
set.seed(505)

for (i in 1:nrow(cenarios)) {
  p <- cenarios$porc[i]
  m <- cenarios$mech[i]
  cat(sprintf("Amputando: Mecanismo %s - Proporção %s...\n", m, p))

  t0 <- Sys.time()

  bancos_amputados <- replicate(1000, rodar_sim(porc = p, metodo = m, dados_base = banco_temp, padrao_mat = matriz_padrao), simplify = FALSE)

  t1 <- Sys.time()
  registro_tempos$tempo_amputacao[i] <- as.numeric(difftime(t1, t0, units = "mins"))

  nome_arq <- sprintf("Resultados/Amputados_%s_%s.rds", m, p)
  saveRDS(bancos_amputados, nome_arq)

  rm(bancos_amputados)
  gc()
}

# Deleção Listwise
cat("\nIniciando Fase 2: Listwise...\n")

cl <- makeCluster(4)
clusterEvalQ(cl, { library(survival) })
clusterExport(cl, "betas_reais")

for (i in 1:nrow(cenarios)) {
  p <- cenarios$porc[i]
  m <- cenarios$mech[i]
  cat(sprintf("Listwise: Mecanismo %s - Proporção %s...\n", m, p))

  nome_arq_amp <- sprintf("Resultados/Amputados_%s_%s.rds", m, p)
  bancos_amputados <- readRDS(nome_arq_amp)

  t0 <- Sys.time()

  resultados_listwise <- pblapply(bancos_amputados, ajustar_listwise_unico, betas_reais_ref = betas_reais, cl = cl)

  t1 <- Sys.time()
  registro_tempos$tempo_listwise[i] <- as.numeric(difftime(t1, t0, units = "mins"))

  nome_arq_lw <- sprintf("Resultados/Listwise_%s_%s.rds", m, p)
  saveRDS(resultados_listwise, nome_arq_lw)

  rm(bancos_amputados, resultados_listwise)
  gc()
}

stopCluster(cl)

# Imputação Múltipla
set.seed(505)
cat("\nIniciando Fase 3: Imputação Múltipla (MICE)...\n")

cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(survival)
  library(mice)
})

for (i in 1:nrow(cenarios)) {
  p <- cenarios$porc[i]
  m <- cenarios$mech[i]
  cat(sprintf("MICE: Mecanismo %s - Proporção %s...\n", m, p))

  nome_arq_amp <- sprintf("Resultados/Amputados_%s_%s.rds", m, p)
  bancos_amputados <- readRDS(nome_arq_amp)

  t0 <- Sys.time()

  resultados_mice <- pblapply(bancos_amputados, ajustar_mice_unico, num_imp = 10, cl = cl)

  t1 <- Sys.time()
  registro_tempos$tempo_mice[i] <- as.numeric(difftime(t1, t0, units = "mins"))

  nome_arq_mice <- sprintf("Resultados/MICE_%s_%s.rds", m, p)
  saveRDS(resultados_mice, nome_arq_mice)

  rm(bancos_amputados, resultados_mice)
  gc()
}

stopCluster(cl)

saveRDS(registro_tempos, "Resultados/Registro_de_Tempos.rds")
print(registro_tempos)
cat("\nProcessamento concluído com sucesso!")