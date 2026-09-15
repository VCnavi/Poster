library(survival)

variaveis <- c("IDADE", "ECGRUP", "GLEASON", "PSA", "CATEATEND", 
               "RECIDIVA", "METASTASE", "QUIMIO", "RADIO", 
               "HORMONIO", "CIRURGIA", "ESCOLARI", "LEI_DIAGTRAT")

modelo_nulo <- survreg(Surv(T_DIAG_O, STATUS) ~ 1, 
                       data = dados_info, 
                       dist = "loglogis")

resultados <- list()

for(var in variaveis) {
  
  formula_uni <- as.formula(paste("Surv(T_DIAG_O, STATUS) ~", var))
  modelo_uni <- survreg(formula_uni, data = dados_info, dist = "loglogis")
  teste_anova <- anova(modelo_nulo, modelo_uni)
  p_valor <- teste_anova$`Pr(>Chi)`[2]
  
  resultados[[var]] <- data.frame(
    Variavel = var,
    P_valor_LRT = p_valor
  )
}

tabela_univariada <- do.call(rbind, resultados)
rownames(tabela_univariada) <- NULL

tabela_univariada$Significativo_5pct <- ifelse(tabela_univariada$P_valor_LRT < 0.05, "Sim", "Não")
tabela_univariada$P_valor_Formatado <- format(tabela_univariada$P_valor_LRT, scientific = FALSE, digits = 4)

tabela_univariada <- tabela_univariada[order(tabela_univariada$P_valor_LRT), ]

modelo_sem_radio <- survreg(
  Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + PSA + CATEATEND +
    RECIDIVA + METASTASE + QUIMIO + HORMONIO + CIRURGIA + ESCOLARI +
    LEI_DIAGTRAT,
  data = dados_info, 
  dist = "loglogis"
)

# A função drop1 para testar a remoção de cada variável
resultado_drop1 <- drop1(modelo_sem_radio, test = "Chisq")

tabela_pvalores <- as.data.frame(resultado_drop1)
tabela_pvalores <- data.frame(
  Variavel = rownames(tabela_pvalores),
  LRT = tabela_pvalores$LRT,
  Pr_Chi = tabela_pvalores$`Pr(>Chi)`
)

tabela_pvalores <- tabela_pvalores[tabela_pvalores$Variavel != "<none>", ]

tabela_pvalores$Significativo <- ifelse(tabela_pvalores$Pr_Chi < 0.05, "Sim", "Não")
tabela_pvalores$P_valor_Formatado <- format(tabela_pvalores$Pr_Chi, scientific = FALSE, digits = 4)

tabela_pvalores <- tabela_pvalores[order(tabela_pvalores$Pr_Chi), ]

print(tabela_pvalores[, c("Variavel", "LRT", "P_valor_Formatado", "Significativo")], row.names = FALSE)

modelo_sem_radio <- survreg(
  Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + PSA + CATEATEND +
    RECIDIVA + METASTASE + QUIMIO + HORMONIO + ESCOLARI +
    LEI_DIAGTRAT,
  data = dados_info, 
  dist = "loglogis"
)

resultado_drop1 <- drop1(modelo_sem_radio, test = "Chisq")

tabela_pvalores <- as.data.frame(resultado_drop1)

# Limpar e formatar a tabela
tabela_pvalores <- data.frame(
  Variavel = rownames(tabela_pvalores),
  LRT = tabela_pvalores$LRT,
  Pr_Chi = tabela_pvalores$`Pr(>Chi)`
)

tabela_pvalores <- tabela_pvalores[tabela_pvalores$Variavel != "<none>", ]

tabela_pvalores$Significativo <- ifelse(tabela_pvalores$Pr_Chi < 0.05, "Sim", "Não")
tabela_pvalores$P_valor_Formatado <- format(tabela_pvalores$Pr_Chi, scientific = FALSE, digits = 4)

tabela_pvalores <- tabela_pvalores[order(tabela_pvalores$Pr_Chi), ]

print(tabela_pvalores[, c("Variavel", "LRT", "P_valor_Formatado", "Significativo")], row.names = FALSE)

modelo_sem_radio <- survreg(
  Surv(T_DIAG_O, STATUS) ~ IDADE + ECGRUP + GLEASON + CATEATEND +
    RECIDIVA + METASTASE + QUIMIO + ESCOLARI +
    LEI_DIAGTRAT,
  data = dados_info, 
  dist = "loglogis"
)

resultado_drop1 <- drop1(modelo_sem_radio, test = "Chisq")

tabela_pvalores <- as.data.frame(resultado_drop1)

tabela_pvalores <- data.frame(
  Variavel = rownames(tabela_pvalores),
  LRT = tabela_pvalores$LRT,
  Pr_Chi = tabela_pvalores$`Pr(>Chi)`
)

tabela_pvalores <- tabela_pvalores[tabela_pvalores$Variavel != "<none>", ]

tabela_pvalores$Significativo <- ifelse(tabela_pvalores$Pr_Chi < 0.05, "Sim", "Não")
tabela_pvalores$P_valor_Formatado <- format(tabela_pvalores$Pr_Chi, scientific = FALSE, digits = 4)

tabela_pvalores <- tabela_pvalores[order(tabela_pvalores$Pr_Chi), ]

print(tabela_pvalores[, c("Variavel", "LRT", "P_valor_Formatado", "Significativo")], row.names = FALSE)