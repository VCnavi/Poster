library(tidyverse)
library(survival)
library(mice)
library(DescTools)
library(broom)
library(parallel)
library(pbapply)

# 1. Preparação dos Valores Reais (Modelo Completo) --------------------
dados <- read.csv("dados_prostata_completos.csv")

modelo_completo <- survreg(
  Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + CATEATEND +
    RECIDIVA + METASTASE + QUIMIO + ESCOLARI + LEI_DIAGTRAT,
  data = dados, 
  dist = "loglogis"
)

tempo_predicao <- 730
betas_reais <- coef(modelo_completo)

# Para LRT e R2 do modelo completo
modelo_nulo_completo <- survreg(Surv(T_DIAG_O, STATUS) ~ 1, data = dados, dist = "loglogis")
loglik_full_real <- modelo_completo$loglik[2]
loglik_null_real <- modelo_nulo_completo$loglik[2]

# Predição média do modelo real em 2 anos
# S(t) para log-logística: 1 / (1 + (t / exp(X*beta))^(1/scale))
escala_real <- modelo_completo$scale
pred_xb_real <- predict(modelo_completo, type = "lp")
surv_2anos_real <- mean(1 / (1 + (tempo_predicao / exp(pred_xb_real))^(1/escala_real)))


# 2. Funções de Extração -----------------------------------------------
# LISTWISE -------------------------------------------------------------
extrair_metricas_listwise <- function(resultado_lw, betas_reais, surv_real, t_pred) {
  
  mod <- resultado_lw$modelo
  
  if (length(mod) == 1 && is.na(mod)) return(NULL) 
  
  # Coeficientes, Wald, P-valor, SE, CI
  resumo <- summary(mod)
  tabela <- as.data.frame(resumo$table)
  colnames(tabela) <- c("Estimativa", "SE", "Wald_Z", "P_valor")
  
  tabela$Variavel <- rownames(tabela)
  tabela <- tabela[!grepl("Log\\(scale\\)", tabela$Variavel), ]
  
  ci <- confint(mod)
  ci <- ci[match(tabela$Variavel, rownames(ci)), , drop = FALSE]
  tabela$CI_Inf <- ci[, 1]
  tabela$CI_Sup <- ci[, 2]
  tabela$AW <- tabela$CI_Sup - tabela$CI_Inf
  
  # Viés e Erro Quadrático 
  tabela$Beta_Real <- betas_reais[tabela$Variavel]
  tabela$Vies <- tabela$Estimativa - tabela$Beta_Real
  tabela$Vies_Perc <- (tabela$Vies / tabela$Beta_Real) * 100
  tabela$Erro_Quad <- (tabela$Estimativa - tabela$Beta_Real)^2
  
  # Cobertura
  tabela$Cobertura <- ifelse(tabela$Beta_Real >= tabela$CI_Inf & tabela$Beta_Real <= tabela$CI_Sup, 1, 0)
  
  # -----------------------------------------------------------
  
  # Teste de Razão de Verossimilhança Global (LRT)
  resposta <- model.response(model.frame(mod))
  mod_nulo <- survreg(resposta ~ 1, dist = "loglogis")
  
  lrt_stat <- 2 * (as.numeric(logLik(mod)) - as.numeric(logLik(mod_nulo)))
  df_diff <- length(coef(mod)) - length(coef(mod_nulo))
  lrt_pval <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
  
  # Predição 2 anos
  escala <- mod$scale
  pred_xb <- predict(mod, type = "lp")
  surv_2anos <- mean(1 / (1 + (t_pred / exp(pred_xb))^(1/escala)))
  diff_pred <- surv_2anos - surv_real
  
  list(
    Tabela_Coefs = tabela,
    Global = data.frame(
      LRT_Pvalor = lrt_pval,
      Surv_2Anos = surv_2anos,
      Diff_Pred_Real = diff_pred
    )
  )
}

# MICE -------------------------------------------------------------------
extrair_metricas_mice <- function(resultado_mice, betas_reais, surv_real, t_pred) {
  
  mids_obj <- resultado_mice$dados_imputados

  mira_completo <- resultado_mice$modelos_brutos 
  mipo_obj <- resultado_mice$modelo_pool
  
  # Coeficientes
  tabela_pool <- summary(mipo_obj, conf.int = TRUE)
  
  colnames(tabela_pool)[colnames(tabela_pool) == "term"] <- "Variavel"
  colnames(tabela_pool)[colnames(tabela_pool) == "estimate"] <- "Estimativa"
  colnames(tabela_pool)[colnames(tabela_pool) == "std.error"] <- "SE"
  colnames(tabela_pool)[colnames(tabela_pool) == "statistic"] <- "Wald_T"
  colnames(tabela_pool)[colnames(tabela_pool) == "p.value"] <- "P_valor"
  colnames(tabela_pool)[colnames(tabela_pool) == "2.5 %"] <- "CI_Inf"
  colnames(tabela_pool)[colnames(tabela_pool) == "97.5 %"] <- "CI_Sup"
  
  tabela_pool$Variavel <- as.character(tabela_pool$Variavel)
  tabela_pool <- tabela_pool[!grepl("Log\\(scale\\)", tabela_pool$Variavel), , drop = FALSE]
  
  # Amplitude do Intervalo
  tabela_pool$AW <- tabela_pool$CI_Sup - tabela_pool$CI_Inf
  
  # Viés e Erro Quadrático
  tabela_pool$Beta_Real <- betas_reais[tabela_pool$Variavel]
  tabela_pool$Vies <- tabela_pool$Estimativa - tabela_pool$Beta_Real
  tabela_pool$Vies_Perc <- (tabela_pool$Vies / tabela_pool$Beta_Real) * 100
  tabela_pool$Erro_Quad <- (tabela_pool$Estimativa - tabela_pool$Beta_Real)^2
  
  # Cobertura
  tabela_pool$Cobertura <- ifelse(tabela_pool$Beta_Real >= tabela_pool$CI_Inf & tabela_pool$Beta_Real <= tabela_pool$CI_Sup, 1, 0)
  
  # -------------------------------------------------------
  
  # Global
  
  mira_nulo <- with(mids_obj, survreg(Surv(T_DIAG_O, STATUS) ~ 1, dist = "loglogis"))
  
  # Wald Multivariado
  teste_global <- tryCatch({
    mice::D1(mira_completo, mira_nulo)
  }, error = function(e) {
    return(NULL)
  })
  
  pval_global <- teste_global$result[, "P(>F)"]
  
  # if (!is.null(teste_global)) {
  #   # Se o D1 funcionou
  #   pval_global <- teste_global$result[, "P(>F)"]
  # } else {
  #   # Se o D1 falhou, calcula o LRT para cada imputação e tira a mediana
  #   lrt_lista <- numeric()
  #   for (mod in mira_completo$analyses) {
  #     resposta <- model.response(model.frame(mod))
  #     mod_nulo_ind <- survreg(resposta ~ 1, dist = "loglogis")
  #     lrt_stat <- 2 * (as.numeric(logLik(mod)) - as.numeric(logLik(mod_nulo_ind)))
  #     df_diff <- length(coef(mod)) - length(coef(mod_nulo_ind))
  #     lrt_p <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
  #     lrt_lista <- c(lrt_lista, lrt_p)
  #   }
  #   pval_global <- median(lrt_lista, na.rm = TRUE)
  # }
  
  # Predições de sobrevivência
  surv_lista <- numeric()
  
  modelos_ajustados <- mira_completo$analyses 
  
  for (mod in modelos_ajustados) {
    escala <- mod$scale
    pred_xb <- predict(mod, type = "lp")
    surv_i <- mean(1 / (1 + (t_pred / exp(pred_xb))^(1 / escala)), na.rm = TRUE)
    surv_lista <- c(surv_lista, surv_i)
  }
  
  surv_2anos_medio <- mean(surv_lista, na.rm = TRUE)
  diff_pred <- surv_2anos_medio - surv_real
  
  list(
    Tabela_Coefs = tabela_pool,
    Global = data.frame(
      Global_Pvalor = pval_global,
      Surv_2Anos = surv_2anos_medio,
      Diff_Pred_Real = diff_pred
    )
  )
}

# Extração -------------------------------------------------------
# Cenários
cenarios <- expand.grid(
  porc = c(0.02),
  mech = c("MAR"),
  stringsAsFactors = FALSE
)

cat("Iniciando extração de métricas...\n")

cl <- makeCluster(4)

clusterEvalQ(cl, {
  library(survival)
  library(mice)
  library(broom)
})

clusterExport(cl, c("extrair_metricas_listwise", 
                    "extrair_metricas_mice", 
                    "betas_reais", 
                    "surv_2anos_real", 
                    "tempo_predicao"))

cat("Iniciando extração de métricas em paralelo...\n")

for (i in 1:nrow(cenarios)) {
  p <- cenarios$porc[i]
  m <- cenarios$mech[i]
  
  cat(sprintf("\nProcessando: Mecanismo %s - Proporção %s\n", m, p))
  
  # --- LISTWISE ---
  arq_lw <- sprintf("Resultados/Listwise_%s_%s.rds", m, p)
  if (file.exists(arq_lw)) {
    dados_lw <- readRDS(arq_lw)
    
    cat("Extraindo métricas (Listwise)...\n")
    metricas_lw <- pblapply(dados_lw, extrair_metricas_listwise, 
                            betas_reais = betas_reais, 
                            surv_real = surv_2anos_real, 
                            t_pred = tempo_predicao,
                            cl = cl)
    
    saveRDS(metricas_lw, sprintf("Resultados/Metricas_Listwise_%s_%s.rds", m, p))
    rm(dados_lw, metricas_lw)
    gc()
  } else {
    warning(sprintf("Arquivo NÃO ENCONTRADO: %s", arq_lw))
  }
  
  # --- MICE ---
  arq_mice <- sprintf("Resultados/MICE_%s_%s.rds", m, p)
  if (file.exists(arq_mice)) {
    dados_mice <- readRDS(arq_mice)
    
    cat("Extraindo métricas (MICE)...\n")
    metricas_mice <- pblapply(dados_mice, extrair_metricas_mice, 
                              betas_reais = betas_reais, 
                              surv_real = surv_2anos_real, 
                              t_pred = tempo_predicao,
                              cl = cl)
    
    saveRDS(metricas_mice, sprintf("Resultados/Metricas_MICE_%s_%s.rds", m, p))
    rm(dados_mice, metricas_mice)
    gc()
  } else {
    warning(sprintf("Arquivo NÃO ENCONTRADO: %s", arq_mice))
  }
}

stopCluster(cl)
cat("\nExtração finalizada. Arquivos salvos na pasta 'Resultados/'.\n")