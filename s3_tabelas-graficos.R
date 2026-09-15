# Compilar resultados -----------------------------------------
# Listwise ----------------------------------------------------
compilar_listwise <- function(cenarios) {
  df_coefs_total <- data.frame()
  df_global_total <- data.frame()
  
  for (i in 1:nrow(cenarios)) {
    p <- cenarios$porc[i]
    m <- cenarios$mech[i]
    arquivo <- sprintf("Resultados/Metricas_Listwise_%s_%s.rds", m, p)
    
    if (file.exists(arquivo)) {
      lista_sim <- readRDS(arquivo)
      
      for (sim in 1:length(lista_sim)) {
        if (is.null(lista_sim[[sim]])) next
        
        # Coeficientes
        coefs <- lista_sim[[sim]]$Tabela_Coefs
        coefs$Mecanismo <- m
        coefs$Proporcao <- p
        coefs$Metodo <- "Listwise"
        coefs$Simulacao <- sim
        
        # Métricas globais
        glob <- lista_sim[[sim]]$Global
        glob$Mecanismo <- m
        glob$Proporcao <- p
        glob$Metodo <- "Listwise"
        glob$Simulacao <- sim
        
        df_coefs_total <- bind_rows(df_coefs_total, coefs)
        df_global_total <- bind_rows(df_global_total, glob)
      }
    }
  }
  return(list(Coefs = df_coefs_total, Global = df_global_total))
}

# MICE --------------------------------------------------------------
compilar_mice <- function(cenarios) {
  df_coefs_total <- data.frame()
  df_global_total <- data.frame()
  
  for (i in 1:nrow(cenarios)) {
    p <- cenarios$porc[i]
    m <- cenarios$mech[i]
    arquivo <- sprintf("Resultados/Metricas_MICE_%s_%s.rds", m, p)
    
    if (file.exists(arquivo)) {
      lista_sim <- readRDS(arquivo)
      
      for (sim in 1:length(lista_sim)) {
        if (is.null(lista_sim[[sim]])) next
        
        # Coeficientes
        coefs <- lista_sim[[sim]]$Tabela_Coefs
        coefs$Mecanismo <- m
        coefs$Proporcao <- p
        coefs$Metodo <- "MICE"
        coefs$Simulacao <- sim
        
        # Métricas globais
        glob <- lista_sim[[sim]]$Global
        glob$Mecanismo <- m
        glob$Proporcao <- p
        glob$Metodo <- "MICE"
        glob$Simulacao <- sim
        
        df_coefs_total <- bind_rows(df_coefs_total, coefs)
        df_global_total <- bind_rows(df_global_total, glob)
      }
    }
  }
  return(list(Coefs = df_coefs_total, Global = df_global_total))
}

# Juntando os dados de Listwise e MICE
cenarios <- expand.grid(porc = c(0.05, 0.1, 0.2, 0.4, 0.6), mech = c("MCAR","MAR","MNAR"))
resultados_lw <- compilar_listwise(cenarios)
resultados_mice <- compilar_mice(cenarios)

df_coefs <- bind_rows(resultados_lw$Coefs, resultados_mice$Coefs)
df_global <- bind_rows(resultados_lw$Global, resultados_mice$Global)

# Tabelas ---------------------------------------------------------
tabela_inferencia <- df_coefs |>
  mutate(
    Vies_Bruto = Estimativa - Beta_Real,
    Nao_Rejeita_Wald = ifelse(P_valor > 0.05, 1, 0)
  ) |>
  group_by(Variavel, Metodo, Mecanismo, Proporcao) |>
  summarise(
    N_Sim = n(),
    Vies_Bruto_Medio = mean(Vies_Bruto, na.rm = TRUE),
    Vies_Rel_Medio = mean(Vies_Perc, na.rm = TRUE),
    Taxa_Nao_Rejeicao = mean(Nao_Rejeita_Wald, na.rm = TRUE) * 100,
    
    EQM = mean(Erro_Quad, na.rm = TRUE),
    
    SD_Vies_Rel = sd(Vies_Perc, na.rm = TRUE),
    SE_Vies_Rel = SD_Vies_Rel / sqrt(N_Sim),
    IC_Vies_Inf = Vies_Rel_Medio - 1.96 * SE_Vies_Rel,
    IC_Vies_Sup = Vies_Rel_Medio + 1.96 * SE_Vies_Rel,
    
    Taxa_Cobertura = mean(Cobertura, na.rm = TRUE),
    SE_Cobertura = sqrt((Taxa_Cobertura * (1 - Taxa_Cobertura)) / N_Sim),
    IC_Cob_Inf = Taxa_Cobertura - 1.96 * SE_Cobertura,
    IC_Cob_Sup = Taxa_Cobertura + 1.96 * SE_Cobertura,
    
    Media_AW = mean(AW, na.rm = TRUE),
    
    .groups = "drop"
  )

# Tabela Global
tabela_predicao <- df_global |>
  group_by(Metodo, Mecanismo, Proporcao) |>
  summarise(
    N_Sim = n(),
    
    # 2 Anos
    Predicao_Media_2A = mean(Surv_2Anos, na.rm = TRUE),
    SE_Pred_2A = sd(Surv_2Anos, na.rm = TRUE) / sqrt(N_Sim),
    IC_Pred_Inf_2A = Predicao_Media_2A - 1.96 * SE_Pred_2A,
    IC_Pred_Sup_2A = Predicao_Media_2A + 1.96 * SE_Pred_2A,
    
    # 4 Anos
    Predicao_Media_4A = mean(Surv_4Anos, na.rm = TRUE),
    SE_Pred_4A = sd(Surv_4Anos, na.rm = TRUE) / sqrt(N_Sim),
    IC_Pred_Inf_4A = Predicao_Media_4A - 1.96 * SE_Pred_4A,
    IC_Pred_Sup_4A = Predicao_Media_4A + 1.96 * SE_Pred_4A,
    
    # 6 Anos
    Predicao_Media_6A = mean(Surv_6Anos, na.rm = TRUE),
    SE_Pred_6A = sd(Surv_6Anos, na.rm = TRUE) / sqrt(N_Sim),
    IC_Pred_Inf_6A = Predicao_Media_6A - 1.96 * SE_Pred_6A,
    IC_Pred_Sup_6A = Predicao_Media_6A + 1.96 * SE_Pred_6A,
    
    # 8 Anos
    Predicao_Media_8A = mean(Surv_8Anos, na.rm = TRUE),
    SE_Pred_8A = sd(Surv_8Anos, na.rm = TRUE) / sqrt(N_Sim),
    IC_Pred_Inf_8A = Predicao_Media_8A - 1.96 * SE_Pred_8A,
    IC_Pred_Sup_8A = Predicao_Media_8A + 1.96 * SE_Pred_8A,
    
    .groups = "drop"
  )

# Gráficos --------------------------------------------------------------
# var_foco seleciona a variável para visualização
var_foco <- "IDADE"
dados_grafico_coefs <- tabela_inferencia |> filter(Variavel == var_foco)
dados_boxplot <- df_coefs |> filter(Variavel == var_foco)

dados <- read.csv("dados_prostata_completos.csv")

modelo_completo <- survreg(
  Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + CATEATEND +
    RECIDIVA + METASTASE + QUIMIO + ESCOLARI + LEI_DIAGTRAT,
  data = dados, 
  dist = "loglogis"
)

# Cálculo dos valores de referência (Modelo Completo)
escala_real <- modelo_completo$scale
pred_xb_real <- predict(modelo_completo, type = "lp")

# Calculando a predição média para cada período (convertido em dias)
surv_2anos_real <- mean(1 / (1 + ((2 * 365) / exp(pred_xb_real))^(1/escala_real)))
surv_4anos_real <- mean(1 / (1 + ((4 * 365) / exp(pred_xb_real))^(1/escala_real)))
surv_6anos_real <- mean(1 / (1 + ((6 * 365) / exp(pred_xb_real))^(1/escala_real)))
surv_8anos_real <- mean(1 / (1 + ((8 * 365) / exp(pred_xb_real))^(1/escala_real)))

# Viés Bruto
g_vies_bruto <- ggplot(dados_grafico_coefs, aes(x = Proporcao, y = Vies_Bruto_Medio, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Mecanismo) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = NULL, 
    y = "Viés Bruto Médio", 
    x = "Proporção de Dados Faltantes"
  ) +
  theme_minimal()

# Viés Relativo
g_vies_relativo <- ggplot(dados_grafico_coefs, aes(x = Proporcao, y = Vies_Rel_Medio, color = Metodo, group = Metodo)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = IC_Vies_Inf, ymax = IC_Vies_Sup), width = 0.02) +
  facet_wrap(~ Mecanismo) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = NULL, y = "Viés Relativo (%)", x = "Proporção de Dados Faltantes") +
  theme_minimal()

# Boxplots das Estimativas
beta_real_val <- unique(dados_boxplot$Beta_Real)[1]

g_boxplot <- ggplot(dados_boxplot, aes(x = as.factor(Proporcao), y = Estimativa, fill = Metodo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + # Omitindo outliers extremos para clareza
  facet_wrap(~ Mecanismo) +
  geom_hline(aes(yintercept = beta_real_val, linetype = "Referência"), color = "red", size = 1) +
  scale_linetype_manual(name = NULL, values = c("Referência" = "dashed"), labels = paste("Referência:", round(beta_real_val, 3))) +
  labs(title = NULL, y = "Estimativa do Coeficiente", x = "Proporção de Dados Faltantes") +
  theme_minimal()

# Taxa de Cobertura
g_cobertura <- ggplot(dados_grafico_coefs, aes(x = Proporcao, y = Taxa_Cobertura, color = Metodo, group = Metodo)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = IC_Cob_Inf, ymax = IC_Cob_Sup), width = 0.02) +
  facet_wrap(~ Mecanismo) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = NULL, y = "Proporção de Cobertura", x = "Proporção de Dados Faltantes") +
  theme_minimal()

# Predição da Sobrevivência em 2 anos 
g_predicao_2A <- ggplot(tabela_predicao, aes(x = Proporcao, y = Predicao_Media_2A, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = IC_Pred_Inf_2A, ymax = IC_Pred_Sup_2A), width = 0.02) +
  facet_wrap(~ Mecanismo) +
  geom_hline(aes(yintercept = surv_2anos_real, linetype = "Referência"), color = "black") + 
  scale_linetype_manual(name = NULL, values = c("Referência" = "dashed"), labels = paste("Referência:", round(surv_2anos_real, 3))) +
  labs(
    title = NULL, 
    y = "Probabilidade Média de Sobrevivência", 
    x = "Proporção de Dados Faltantes"
  ) +
  theme_minimal()

# Predição da Sobrevivência em 4 anos
g_predicao_4A <- ggplot(tabela_predicao, aes(x = Proporcao, y = Predicao_Media_4A, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = IC_Pred_Inf_4A, ymax = IC_Pred_Sup_4A), width = 0.02) +
  facet_wrap(~ Mecanismo) +
  geom_hline(aes(yintercept = surv_4anos_real, linetype = "Referência"), color = "black") + 
  scale_linetype_manual(name = NULL, values = c("Referência" = "dashed"), labels = paste("Referência:", round(surv_4anos_real, 3))) +
  labs(
    title = NULL, 
    y = "Probabilidade Média de Sobrevivência", 
    x = "Proporção de Dados Faltantes"
  ) +
  theme_minimal()

# Predição da Sobrevivência em 6 anos
g_predicao_6A <- ggplot(tabela_predicao, aes(x = Proporcao, y = Predicao_Media_6A, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = IC_Pred_Inf_6A, ymax = IC_Pred_Sup_6A), width = 0.02) +
  facet_wrap(~ Mecanismo) +
  geom_hline(aes(yintercept = surv_6anos_real, linetype = "Referência"), color = "black") + 
  scale_linetype_manual(name = NULL, values = c("Referência" = "dashed"), labels = paste("Referência:", round(surv_6anos_real, 3))) +
  labs(
    title = NULL, 
    y = "Probabilidade Média de Sobrevivência", 
    x = "Proporção de Dados Faltantes"
  ) +
  theme_minimal()

# Predição da Sobrevivência em 8 anos
g_predicao_8A <- ggplot(tabela_predicao, aes(x = Proporcao, y = Predicao_Media_8A, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = IC_Pred_Inf_8A, ymax = IC_Pred_Sup_8A), width = 0.02) +
  facet_wrap(~ Mecanismo) +
  geom_hline(aes(yintercept = surv_8anos_real, linetype = "Referência"), color = "black") + 
  scale_linetype_manual(name = NULL, values = c("Referência" = "dashed"), labels = paste("Referência:", round(surv_8anos_real, 3))) +
  labs(
    title = NULL, 
    y = "Probabilidade Média de Sobrevivência", 
    x = "Proporção de Dados Faltantes"
  ) +
  theme_minimal()

# Porcentagem de não rejeição de Wald
g_wald <- ggplot(dados_grafico_coefs, aes(x = Proporcao, y = Taxa_Nao_Rejeicao, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Mecanismo) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = NULL, 
    y = "Não Rejeição de H0 (%)", 
    x = "Proporção de Dados Faltantes"
  ) +
  theme_minimal()

# Erro Quadrático Médio
g_eqm <- ggplot(dados_grafico_coefs, aes(x = Proporcao, y = EQM, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Mecanismo) +
  labs(
    title = NULL, 
    y = "EQM", 
    x = "Proporção de Dados Faltantes"
  ) +
  theme_minimal()

# Amplitude
g_amplitude <- ggplot(dados_grafico_coefs, aes(x = Proporcao, y = Media_AW, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Mecanismo) +
  labs(title = NULL, y = "Amplitude", x = "Proporção de Dados Faltantes") +
  theme_minimal()

# Exibir os gráficos
g_vies_bruto
g_vies_relativo
g_boxplot
g_cobertura
g_predicao_2A
g_predicao_4A
g_predicao_6A
g_predicao_8A
g_wald
g_amplitude
g_eqm
