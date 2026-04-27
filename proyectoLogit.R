# ==============================================================================
# PROYECTO DE REGRESIÓN LOGÍSTICA
# DANIELA RENÉE Y MARIO SÁNCHEZ
# ECONOMETRÍA 1
# ==============================================================================
#Ruta Dani Abajo
##setwd("/Users/danielarenee/Desktop/Water_Quality_Logit")
#Ruta Mario Abajo
setwd("C:/Users/msfcl/OneDrive/Escritorio/Water_Quality_Logit")

library(car)      
library(corrplot)
library(dplyr)


PATH_TRAIN   <- "base_train.csv"
PATH_TEST    <- "base_test.csv"
PATH_OUT     <- "base_train_final.csv"
PATH_OUT_TS  <- "base_test_final.csv"
PATH_HEATMAP <- "heatmap_correlaciones_R.png"
PATH_HEATMAP2 <- "heatmap_correlaciones2_R.png"

UMBRAL_VIF   <- 10   # VIF > 10 se considera multicolinealidad crítica

# variables que no entran al modelo
METADATOS <- c("CLAVE SITIO", "CLAVE DE MONITOREO", "NOMBRE DEL SITIO",
               "FECHA REALIZACIÓN", "Año",
               "DQO_ESTADO", "SST_ESTADO", "EC_ESTADO", "SEMAFORO",
               "TIPO CUERPO DE AGUA")

ELIMINAR_MANUAL <- c("SDT", "OD_.") 

# CARGA DE DATOS
cat("1: CARGA DE DATOS\n")

train <- read.csv(PATH_TRAIN, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")
test  <- read.csv(PATH_TEST,  check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")

# renombrar OD_% y OD mg/Lpara R
names(train)[names(train) == "OD_%"] <- "OD_pct"
names(test)[names(test)  == "OD_%"]  <- "OD_pct"
names(train)[names(train) == "OD_mg/L"] <- "OD_mgL"
names(test)[names(test)   == "OD_mg/L"] <- "OD_mgL"

ELIMINAR_MANUAL <- c("SDT", "OD_pct")

cat(sprintf("  Train: %d filas x %d columnas\n", nrow(train), ncol(train)))
cat(sprintf("  Test:  %d filas x %d columnas\n", nrow(test),  ncol(test)))
cat(sprintf("  Prevalencia y=1 en train: %.3f\n", mean(train$y)))
cat(sprintf("  Prevalencia y=1 en test:  %.3f\n", mean(test$y)))

# IDENTIFICAR CANDIDATAS 
cat("2: POOL INICIAL DE CANDIDATAS\n")

# candidatas = todas las numéricas menos metadatos, y, y las que eliminamos a mano
cols_excluir <- c(METADATOS, "y", ELIMINAR_MANUAL)
candidatas   <- setdiff(names(train), cols_excluir)
# solo numéricas
candidatas   <- candidatas[sapply(train[candidatas], is.numeric)]

cat(sprintf("  Candidatas iniciales: %d\n", length(candidatas)))
cat(sprintf("  Eliminadas manualmente: %s\n",
            paste(ELIMINAR_MANUAL, collapse = ", ")))

# MATRIZ DE CORRELACIONES DE SPEARMAN 
cat("3: MATRIZ DE CORRELACIONES\n")

corr <- cor(train[, candidatas], method = "spearman")

# pares con |r| > 0.7
cat("  Pares con |correlación Spearman| > 0.7:\n")
for (i in 1:(ncol(corr) - 1)) {
  for (j in (i + 1):ncol(corr)) {
    if (abs(corr[i, j]) > 0.7) {
      cat(sprintf("    %-20s <-> %-20s  r = %+.3f\n",
                  colnames(corr)[i], colnames(corr)[j], corr[i, j]))
    }
  }
}

# HEATMAP
cat("4: HEATMAP DE CORRELACIONES\n")

png(PATH_HEATMAP, width = 1400, height = 1300, res = 120)
corrplot(corr,
         method    = "color",
         type      = "full",
         order     = "original",
         tl.col    = "black",
         tl.srt    = 90,
         tl.cex    = 0.7,
         addCoef.col = "black",
         number.cex  = 0.45,
         col       = colorRampPalette(c("#016FB9", "white", "#B6244F"))(200),
         title     = "Matriz de correlaciones de Spearman (train, R)",
         mar       = c(0, 0, 2, 0))
dev.off()
cat(sprintf("  ✔ Heatmap guardado en: %s\n", PATH_HEATMAP))

# 5: VIF ITERATIVO
#   - calcular VIF de todas las variables
#   - identificar todas las que tengan VIF > umbral
#   - entre esas, eliminar la que menos correlacione con y (Spearman)
#   - recalcular, repetir hasta que todas tengan VIF <= umbral

cat("5: VIF ITERATIVO\n")

# funciones auxiliares
calcular_vif <- function(vars, data) {
  formula <- as.formula(paste("y ~", paste(paste0("`", vars, "`"),
                                           collapse = " + ")))
  modelo  <- lm(formula, data = data)
  vif(modelo)
}

# correlación de Spearman de cada candidata con y
corr_con_y <- sapply(candidatas, function(v) {
  abs(cor(train[[v]], train$y, method = "spearman"))
})

vars_actuales <- candidatas
iter <- 0

repeat {
  iter <- iter + 1
  
  vifs <- calcular_vif(vars_actuales, train)
  max_vif <- max(vifs)
  
  cat(sprintf("\n  Iteración %d — %d variables | peor VIF: %.2f\n",
              iter, length(vars_actuales), max_vif))
  
  if (max_vif <= UMBRAL_VIF) {
    cat("  ✔ Todas las variables con VIF <= ", UMBRAL_VIF, "\n", sep = "")
    break
  }
  
  # identificar candidatas a eliminar: todas las que tengan VIF > umbral
  culpables <- names(vifs[vifs > UMBRAL_VIF])
  
  # entre las culpables, elegir la de menor correlación con y
  corrs_culpables <- corr_con_y[culpables]
  victima <- names(which.min(corrs_culpables))
  
  cat(sprintf("    Culpables (VIF > %d): %d variables\n",
              UMBRAL_VIF, length(culpables)))
  for (v in culpables) {
    marca <- if (v == victima) " <-- eliminada" else ""
    cat(sprintf("      %-22s VIF=%6.2f  |r_Spearman(y)|=%.3f%s\n",
                v, vifs[v], corr_con_y[v], marca))
  }
  
  vars_actuales <- setdiff(vars_actuales, victima)
}

cat("VIF FINAL DE LAS VARIABLES SUPERVIVIENTES\n")
vifs_finales <- calcular_vif(vars_actuales, train)
vifs_finales_ord <- sort(vifs_finales, decreasing = TRUE)
for (v in names(vifs_finales_ord)) {
  cat(sprintf("  %-22s VIF = %6.3f   |r_Spearman(y)| = %.3f\n",
              v, vifs_finales_ord[v], corr_con_y[v]))
}

cat(sprintf("\n  Variables iniciales: %d\n", length(candidatas)))
cat(sprintf("  Variables finales:   %d\n", length(vars_actuales)))
cat(sprintf("  Eliminadas por VIF:  %d\n",
            length(candidatas) - length(vars_actuales)))

# EXPORTAR POOL
# Conservamos metadatos, las candidatas sobrevivientes, y y
cols_finales <- c(intersect(METADATOS, names(train)), vars_actuales, "y")

train_final <- train[, cols_finales]
test_final  <- test[,  cols_finales]

write.csv(train_final, PATH_OUT,    row.names = FALSE, fileEncoding = "UTF-8")
write.csv(test_final,  PATH_OUT_TS, row.names = FALSE, fileEncoding = "UTF-8")

cat(sprintf("  ✔ Train final: %s  (%d x %d)\n",
            PATH_OUT,    nrow(train_final), ncol(train_final)))
cat(sprintf("  ✔ Test final:  %s  (%d x %d)\n",
            PATH_OUT_TS, nrow(test_final),  ncol(test_final)))
cat(sprintf("\n  Pool final de regresoras: %s\n",
            paste(vars_actuales, collapse = ", ")))


# HEATMAP 2
cat("4: HEATMAP DE CORRELACIONES 2\n")
corr_final <- cor(train[, vars_actuales], method = "spearman")

png(PATH_HEATMAP2, width = 1400, height = 1300, res = 120)
corrplot(corr_final,
         method    = "color",
         type      = "full",
         order     = "original",
         tl.col    = "black",
         tl.srt    = 90,
         tl.cex    = 0.7,
         addCoef.col = "black",
         number.cex  = 0.45,
         col       = colorRampPalette(c("#016FB9", "white", "#B6244F"))(200),
         title     = "Matriz de correlaciones de Spearman Post-VIF (train, R)",
         mar       = c(0, 0, 2, 0))
dev.off()

cat(sprintf("  ✔ Heatmap final guardado en: %s\n", PATH_HEATMAP2))


# AJUSTE DE LOS 4 MODELOS CANDIDATOS

# modelo base (con las 28)
formula_base <- as.formula(paste("y ~", paste(paste0("`", vars_actuales, "`"),
                                            collapse = " + ")))
modelo_base <- glm(formula_base, data = train, family = binomial(link = "logit"))
cat(sprintf("  Variables: %d\n", length(vars_actuales)))
cat(sprintf("  AIC: %.2f   BIC: %.2f\n", AIC(modelo_base), BIC(modelo_base)))


# Proponemos 4 modelos de regresión logística:
#   Modelo 1 : 8 variables basado en la teoría 
#   Modelo 2 (stepwise AIC): selección bidireccional a partir del modelo amplio
#   Modelo 3 (stepwise BIC): selección bidireccional a partir del modelo amplio

cat("6: AJUSTE DE LOS 3 MODELOS CANDIDATOS\n")

# ------------------------------------------------------------------------------
# MODELO 1: Basado en la teoría (8 variables)
# ------------------------------------------------------------------------------
cat("\n--- Modelo 1: Basado en teoría ---\n")

vars_m1 <- c("TEMP_AGUA", "SST", "COLOR_VER",
             "pH_CAMPO", "CONDUC_CAMPO", "OD_mgL", "NI_TOT", "E_COLI")

formula_m1 <- as.formula(paste("y ~", paste(paste0("`", vars_m1, "`"),
                                            collapse = " + ")))
modelo_1 <- glm(formula_m1, data = train, family = binomial(link = "logit"))
cat(sprintf("  Variables: %d\n", length(vars_m1)))
cat(sprintf("  AIC: %.2f   BIC: %.2f\n", AIC(modelo_1), BIC(modelo_1)))

# ------------------------------------------------------------------------------
# MODELO 2: Stepwise por AIC (bidireccional)
# ------------------------------------------------------------------------------
# step() arranca desde modelo_base y en cada paso considera agregar o quitar
cat("\n--- Modelo 2: Stepwise AIC (bidireccional) ---\n")

modelo_2 <- step(modelo_base, direction = "both", trace = 0, k = 2)
vars_m2  <- names(coef(modelo_2))[-1]   # quitar el intercepto
# limpiar backticks del nombre de cada coeficiente
vars_m2  <- gsub("`", "", vars_m2)
cat(sprintf("  Variables seleccionadas: %d\n", length(vars_m2)))
cat(sprintf("  AIC: %.2f   BIC: %.2f\n", AIC(modelo_2), BIC(modelo_2)))
cat("  Variables:", paste(vars_m2, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# MODELO 3: Stepwise por BIC (bidireccional)
# ------------------------------------------------------------------------------
# BIC penaliza más la complejidad que AIC
cat("\n--- Modelo 3: Stepwise BIC (bidireccional) ---\n")

n_train  <- nrow(train)
modelo_3 <- step(modelo_base, direction = "both", trace = 0, k = log(n_train))
vars_m3  <- names(coef(modelo_3))[-1]
vars_m3  <- gsub("`", "", vars_m3)
cat(sprintf("  Variables seleccionadas: %d\n", length(vars_m3)))
cat(sprintf("  AIC: %.2f   BIC: %.2f\n", AIC(modelo_3), BIC(modelo_3)))
cat("  Variables:", paste(vars_m3, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# TABLA RESUMEN COMPARATIVA
# ------------------------------------------------------------------------------
cat("\n--- Tabla resumen de los 3 modelos ---\n")
resumen <- data.frame(
  Modelo   = c("M1 Teoría", "M2 Stepwise AIC", "M3 Stepwise BIC"),
  n_vars   = c(length(vars_m1), length(vars_actuales),
               length(vars_m3)),
  AIC      = c(AIC(modelo_1), AIC(modelo_2), AIC(modelo_3)),
  BIC      = c(BIC(modelo_1), BIC(modelo_2), BIC(modelo_3)),
  deviance = c(modelo_1$deviance, modelo_2$deviance,
               modelo_3$deviance)
)
resumen[, c("AIC", "BIC", "deviance")] <- round(resumen[, c("AIC", "BIC", "deviance")], 2)
print(resumen, row.names = FALSE)


# ==============================================================================
# INFERENCIA DE WALD - PRUEBA DE SIGNIFICANCIA INDIVIDUAL

# Para cada coeficiente del modelo se prueba:
#   H0: beta_j = 0    vs    H1: beta_j != 0
#   Z_0 = beta_hat_j / se(beta_hat_j)
# donde se viene de la diagonal de Var(beta) = -G(beta)^(-1) (inversa negativa
# de la hessiana). Estos son los valores que ya reporta summary(glm).
#
# Bajo H0 y muestras grandes, Z_0 ~ N(0,1), por lo que:
#   p-valor = 2 * (1 - Phi(|Z_0|))
#
# Criterio de decisión con alpha = 0.05:
#   Rechazar H0 si p-valor < 0.05
#
# Códigos de significancia:
#   ***  p < 0.001    evidencia muy fuerte
#   **   p < 0.01     evidencia fuerte
#   *    p < 0.05     evidencia moderada
#   ns   p >= 0.05    no significativo
# ==============================================================================

cat("\n==============================================================\n")
cat("7: INFERENCIA DE WALD POR COEFICIENTE\n")
cat("==============================================================\n")

# función para construir la tabla de Wald de cualquier glm
tabla_wald <- function(modelo, nombre_modelo, alpha = 0.05) {
  
  coefs <- summary(modelo)$coefficients
  # columnas que entrega R: Estimate, Std. Error, z value, Pr(>|z|)
  
  tabla <- data.frame(
    Variable     = rownames(coefs),
    Beta         = round(coefs[, "Estimate"],   4),
    SE           = round(coefs[, "Std. Error"], 4),
    Z            = round(coefs[, "z value"],    3),
    p_valor      = coefs[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  # decisión al alpha = 0.05
  tabla$Decision <- ifelse(tabla$p_valor < alpha,
                           "Rechazar H0", "No rechazar H0")
  
  # códigos de significancia
  tabla$Sig <- cut(tabla$p_valor,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns"),
                   right  = FALSE)
  
  # formato del p-valor para impresión
  tabla$p_valor <- format.pval(tabla$p_valor, digits = 4, eps = 1e-4)
  
  # limpiar backticks de los nombres
  tabla$Variable <- gsub("`", "", tabla$Variable)
  
  cat(sprintf("\n--- %s ---\n", nombre_modelo))
  cat(sprintf("  alpha = %.2f   |   z_critico = %.3f\n",
              alpha, qnorm(1 - alpha/2)))
  print(tabla, row.names = FALSE)
  
  # resumen
  n_rechaza  <- sum(tabla$Decision == "Rechazar H0") - 
    as.integer("(Intercept)" %in% tabla$Variable &
                 tabla$Decision[tabla$Variable == "(Intercept)"] == "Rechazar H0")
  # corrección: contar excluyendo intercepto
  sin_intercepto <- tabla[tabla$Variable != "(Intercept)", ]
  n_rechaza  <- sum(sin_intercepto$Decision == "Rechazar H0")
  n_total    <- nrow(sin_intercepto)
  
  cat(sprintf("\n  Variables significativas (sin contar intercepto): %d de %d (%.1f%%)\n",
              n_rechaza, n_total, 100 * n_rechaza / n_total))
  
  invisible(tabla)
}

# Aplicamos a los 3 modelos
tabla_m1 <- tabla_wald(modelo_1, "Modelo 1: Teórico")
tabla_m2 <- tabla_wald(modelo_2, "Modelo 2: Stepwise AIC")
tabla_m3 <- tabla_wald(modelo_3, "Modelo 3: Stepwise BIC")


#PRUEBA DE RAZÓN DE VEROSIMILITUDES (LRT) PARA VARIABLES CON SEPARACIÓN
# ==============================================================================
# La prueba de Wald NO es confiable cuando una variable presenta separación
# cuasi-completa: cuando predice y casi perfectamente, el estimador de
# máxima verosimilitud beta_hat tiende a +-Inf, el error estándar también
# explota, y el cociente Z = beta/se -> 0/Inf arroja un p-valor cercano a 1
# que NO refleja la verdadera contribución de la variable al modelo.
#
# Detección de separación (criterio mixto):
#
#   (a) Variables BINARIAS (0/1):
#       Se inspecciona la tabla cruzada 2x2 contra y. Si alguna celda es
#       cero, hay separación cuasi-completa por definición. Este criterio
#       es exacto, no heurístico, y no depende de umbrales arbitrarios.
#
#   (b) Variables CONTINUAS:
#       Se usa SE > 30 como bandera (Allison 2008, "Convergence Failures
#       in Logistic Regression", SAS Global Forum, paper 360). En la
#       práctica las continuas casi nunca presentan separación, pero
#       conviene tener el diagnóstico.
#
# Solución estándar (Hosmer & Lemeshow 2013, cap. 4): usar la prueba de
# razón de verosimilitudes (LRT), que es el "método de desviación
# condicional" de las notas de clase. Para H0: beta_j = 0:
#
#   G2 = lambda(beta_j | beta_{-j}) = lambda(beta_{-j}) - lambda(beta)
#      = -2 * [logL(modelo_reducido) - logL(modelo_completo)]
#
# Bajo H0, G2 ~ chi^2 con 1 grado de libertad. Rechazamos H0 si p < alpha.
#
# A diferencia de Wald, LRT NO requiere el error estándar de beta_j, por
# lo que es válido aun bajo separación. Además es la prueba óptima de
# máxima verosimilitud (Wilks 1938).
# ==============================================================================

cat("\n==============================================================\n")
cat("7B: LRT PARA VARIABLES CON SEPARACIÓN CUASI-COMPLETA\n")
cat("==============================================================\n")

# ------------------------------------------------------------------------------
# Detector de separación
# ------------------------------------------------------------------------------
# Devuelve los nombres de las variables del modelo que presentan separación.
detectar_separacion <- function(modelo, datos, umbral_se_continua = 30) {
  
  # variables del modelo (sin intercepto)
  vars_modelo <- attr(terms(modelo), "term.labels")
  vars_modelo <- gsub("`", "", vars_modelo)
  
  con_separacion <- character(0)
  motivo         <- character(0)
  
  coefs <- summary(modelo)$coefficients
  
  for (v in vars_modelo) {
    valores <- datos[[v]]
    
    # caso binaria: usar tabla cruzada
    if (all(valores %in% c(0, 1))) {
      tabla <- table(valores, datos$y)
      if (any(tabla == 0)) {
        con_separacion <- c(con_separacion, v)
        motivo <- c(motivo, "celda 0 en tabla cruzada")
      }
      next
    }
    
    # caso continua: revisar SE
    fila <- coefs[grepl(v, rownames(coefs), fixed = TRUE), , drop = FALSE]
    if (nrow(fila) > 0) {
      se <- fila[1, "Std. Error"]
      if (se > umbral_se_continua) {
        con_separacion <- c(con_separacion, v)
        motivo <- c(motivo, sprintf("SE = %.1f > 30 (Allison 2008)", se))
      }
    }
  }
  
  if (length(con_separacion) == 0) {
    return(NULL)
  }
  data.frame(Variable = con_separacion, Motivo = motivo,
             stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------------
# LRT individual para una lista de variables
# ------------------------------------------------------------------------------
lrt_separacion <- function(modelo, datos_train, nombre_modelo, alpha = 0.05) {
  
  cat(sprintf("\n--- %s ---\n", nombre_modelo))
  
  diag <- detectar_separacion(modelo, datos_train)
  
  if (is.null(diag)) {
    cat("  Ninguna variable presenta separación.\n")
    cat("  Las pruebas de Wald son confiables para todas las variables.\n")
    return(invisible(NULL))
  }
  
  cat("  Variables con separación detectadas:\n")
  for (i in seq_len(nrow(diag))) {
    cat(sprintf("    - %-15s (%s)\n", diag$Variable[i], diag$Motivo[i]))
  }
  cat("  Re-evaluando significancia con LRT (chi-cuadrada, 1 gl):\n\n")
  
  resultados <- data.frame(
    Variable = character(), G2 = numeric(), df = integer(),
    p_valor  = numeric(), Decision = character(),
    stringsAsFactors = FALSE
  )
  
  for (var in diag$Variable) {
    formula_completa <- formula(modelo)
    formula_reducida <- update(formula_completa,
                               as.formula(paste(". ~ . -", paste0("`", var, "`"))))
    modelo_reducido  <- glm(formula_reducida, data = datos_train,
                            family = binomial(link = "logit"))
    
    test <- anova(modelo_reducido, modelo, test = "Chisq")
    G2  <- test$Deviance[2]
    df  <- test$Df[2]
    pv  <- test$`Pr(>Chi)`[2]
    
    resultados <- rbind(resultados, data.frame(
      Variable = var,
      G2       = round(G2, 3),
      df       = df,
      p_valor  = pv,
      Decision = ifelse(pv < alpha, "Rechazar H0", "No rechazar H0"),
      stringsAsFactors = FALSE
    ))
  }
  
  resultados$p_valor <- format.pval(resultados$p_valor, digits = 4, eps = 1e-4)
  print(resultados, row.names = FALSE)
  invisible(resultados)
}

# ------------------------------------------------------------------------------
# Aplicar a los 3 modelos
# ------------------------------------------------------------------------------
lrt_m1 <- lrt_separacion(modelo_1, train, "Modelo 1: Teórico")
lrt_m2 <- lrt_separacion(modelo_2, train, "Modelo 2: Stepwise AIC")
lrt_m3 <- lrt_separacion(modelo_3, train, "Modelo 3: Stepwise BIC")



# ==============================================================================
#  VERSIONES SIN COT_NOM DE M2 Y M3
# ==============================================================================
# COT_NOM presenta separación cuasi-completa contra y porque el Carbono
# Orgánico Total (COT) y la Demanda Química de Oxígeno (DQO, que define y)
# son medidas alternativas del mismo fenómeno físico-químico: la carga
# orgánica del agua. La literatura de química ambiental documenta que
# DQO ≈ 3·COT en muestras típicas (Stumm & Morgan, "Aquatic Chemistry").
#
# Por esta razón, conservar COT_NOM en los modelos finales genera una
# predicción casi tautológica: "el agua incumple la norma en COT (carga
# orgánica) -> el agua incumple la norma en DQO (carga orgánica)". El
# modelo aprende algo trivial, no nuevo.
#
# Construimos versiones "sin COT" de M2 y M3 para evaluar cuánto poder
# predictivo aporta COT_NOM más allá de las demás regresoras y para
# tener modelos finales más informativos sobre las otras dimensiones de
# la calidad del agua (nutrientes, físicos, microbiológicos).
# ==============================================================================

cat("\n==============================================================\n")
cat("6B: VERSIONES SIN COT_NOM\n")
cat("==============================================================\n")

# ------------------------------------------------------------------------------
# MODELO 2 SIN COT: Stepwise AIC sin COT_NOM
# ------------------------------------------------------------------------------
cat("\n--- Modelo 2 sin COT: Stepwise AIC ---\n")

vars_m2_sinCOT <- setdiff(vars_m2, "COT_NOM")
formula_m2_sc  <- as.formula(paste("y ~", paste(paste0("`", vars_m2_sinCOT, "`"),
                                                collapse = " + ")))
modelo_2_sinCOT <- glm(formula_m2_sc, data = train,
                       family = binomial(link = "logit"))
cat(sprintf("  Variables: %d\n", length(vars_m2_sinCOT)))
cat(sprintf("  AIC: %.2f   BIC: %.2f\n",
            AIC(modelo_2_sinCOT), BIC(modelo_2_sinCOT)))
cat("  Variables:", paste(vars_m2_sinCOT, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# MODELO 3 SIN COT: Stepwise BIC sin COT_NOM
# ------------------------------------------------------------------------------
cat("\n--- Modelo 3 sin COT: Stepwise BIC ---\n")

vars_m3_sinCOT <- setdiff(vars_m3, "COT_NOM")
formula_m3_sc  <- as.formula(paste("y ~", paste(paste0("`", vars_m3_sinCOT, "`"),
                                                collapse = " + ")))
modelo_3_sinCOT <- glm(formula_m3_sc, data = train,
                       family = binomial(link = "logit"))
cat(sprintf("  Variables: %d\n", length(vars_m3_sinCOT)))
cat(sprintf("  AIC: %.2f   BIC: %.2f\n",
            AIC(modelo_3_sinCOT), BIC(modelo_3_sinCOT)))
cat("  Variables:", paste(vars_m3_sinCOT, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# COMPARACIÓN: ¿cuánto se pierde al sacar COT_NOM?
# ------------------------------------------------------------------------------
# Si AIC sube poco al sacar COT_NOM, su valor predictivo es marginal y se
# justifica la decisión metodológica de eliminarla. Si sube mucho, es señal
# de que aporta información sustantiva más allá de la redundancia con DQO.
cat("\n--- Comparación con vs sin COT_NOM ---\n")
comp <- data.frame(
  Modelo       = c("M2 con COT", "M2 sin COT",
                   "M3 con COT", "M3 sin COT"),
  n_vars       = c(length(vars_m2),        length(vars_m2_sinCOT),
                   length(vars_m3),        length(vars_m3_sinCOT)),
  AIC          = c(AIC(modelo_2),          AIC(modelo_2_sinCOT),
                   AIC(modelo_3),          AIC(modelo_3_sinCOT)),
  BIC          = c(BIC(modelo_2),          BIC(modelo_2_sinCOT),
                   BIC(modelo_3),          BIC(modelo_3_sinCOT)),
  deviance     = c(modelo_2$deviance,      modelo_2_sinCOT$deviance,
                   modelo_3$deviance,      modelo_3_sinCOT$deviance)
)
comp[, c("AIC","BIC","deviance")] <- round(comp[, c("AIC","BIC","deviance")], 2)
print(comp, row.names = FALSE)

# Cuantificamos cuánto pierde cada modelo al quitar COT_NOM
cat("\n  Penalización por quitar COT_NOM (en AIC):\n")
cat(sprintf("    M2 -> M2 sin COT:  delta AIC = %+.2f\n",
            AIC(modelo_2_sinCOT) - AIC(modelo_2)))
cat(sprintf("    M3 -> M3 sin COT:  delta AIC = %+.2f\n",
            AIC(modelo_3_sinCOT) - AIC(modelo_3)))

# ==============================================================================
# 8: PRUEBAS DE DESVIACIÓN 
# ==============================================================================
#
#
# PRUEBA 1 - Desviación del modelo (vs modelo saturado teórico):
#   lambda(beta) = 2*ln(modsat) - 2*ln(L(beta_hat))
#   Bajo H0: el modelo es correcto, lambda(beta) ~ chi^2_{n-p}
#   Criterio: NO RECHAZAR si lambda(beta) <= chi^2_{alpha, n-p}
#   En R, lambda(beta) = modelo$deviance (la "Residual deviance")
#
# PRUEBA 2 - Subconjuntos de parámetros:
#   H0: beta_2 = 0  vs  H1: beta_2 != 0
#   lambda(beta_2 | beta_1) = lambda(beta_1) - lambda(beta)
#   Bajo H0, distribuye chi^2_r donde r = nro de parámetros omitidos
#   Criterio: RECHAZAR H0 si lambda(beta_2 | beta_1) >= chi^2_{alpha, r}
#
# El "modelo completo" beta para la Prueba 2 es el modelo base post-VIF
# excluyendo COT_NOM (27 variables), por consistencia con los candidatos
# limpios. Cada candidato es el "modelo reducido" beta_1.
# ==============================================================================

cat("\n==============================================================\n")
cat("8: PRUEBAS DE DESVIACIÓN (notas de clase)\n")
cat("==============================================================\n")

ALPHA <- 0.05

# ------------------------------------------------------------------------------
# Modelo base sin COT_NOM (27 variables) - referencia para Prueba 2
# ------------------------------------------------------------------------------
vars_base_sinCOT <- setdiff(vars_actuales, "COT_NOM")
formula_base_sc  <- as.formula(paste("y ~", paste(paste0("`", vars_base_sinCOT, "`"),
                                                  collapse = " + ")))
modelo_base_sinCOT <- glm(formula_base_sc, data = train,
                          family = binomial(link = "logit"))
cat(sprintf("\nModelo base sin COT (referencia): %d variables\n",
            length(vars_base_sinCOT)))
cat(sprintf("  Deviance: %.2f   AIC: %.2f\n",
            modelo_base_sinCOT$deviance, AIC(modelo_base_sinCOT)))

# ------------------------------------------------------------------------------
# PRUEBA 1: Desviación del modelo (vs saturado teórico)
# ------------------------------------------------------------------------------
cat("\n--------------------------------------------------------------\n")
cat("PRUEBA 1: Desviación del modelo lambda(beta)\n")
cat("  H0: el modelo es adecuado (vs saturado teórico)\n")
cat("  Estadístico: lambda(beta) = modelo$deviance\n")
cat("  Distribución bajo H0: chi^2 con n-p grados de libertad\n")
cat("  Criterio: NO RECHAZAR H0 si lambda <= chi^2_{alpha, n-p}\n")
cat("--------------------------------------------------------------\n")

prueba1 <- function(modelo, nombre, n, alpha = ALPHA) {
  lambda  <- modelo$deviance
  p       <- length(coef(modelo))      # parámetros (incluye intercepto)
  gl      <- n - p
  critico <- qchisq(1 - alpha, df = gl)
  pval    <- 1 - pchisq(lambda, df = gl)
  decision <- if (lambda <= critico) "No rechazar H0 (modelo adecuado)"
  else                   "Rechazar H0 (modelo NO adecuado)"
  
  cat(sprintf("\n  %s\n", nombre))
  cat(sprintf("    lambda(beta)        = %.3f\n", lambda))
  cat(sprintf("    gl (n - p)          = %d - %d = %d\n", n, p, gl))
  cat(sprintf("    chi^2 critico (a=%.2f) = %.3f\n", alpha, critico))
  cat(sprintf("    p-valor             = %s\n",
              format.pval(pval, digits = 4, eps = 1e-4)))
  cat(sprintf("    Decision            = %s\n", decision))
  
  data.frame(Modelo = nombre, lambda = round(lambda, 3),
             gl = gl, critico = round(critico, 3),
             p_valor = pval, Decision = decision,
             stringsAsFactors = FALSE)
}

n_train <- nrow(train)
p1_m1  <- prueba1(modelo_1,        "Modelo 1: Teórico",        n_train)
p1_m2s <- prueba1(modelo_2_sinCOT, "Modelo 2 sin COT (AIC)",   n_train)
p1_m3s <- prueba1(modelo_3_sinCOT, "Modelo 3 sin COT (BIC)",   n_train)

# ------------------------------------------------------------------------------
# PRUEBA 2: Subconjuntos de parámetros (vs modelo base sin COT)
# ------------------------------------------------------------------------------
cat("\n\n--------------------------------------------------------------\n")
cat("PRUEBA 2: Subconjuntos de parámetros lambda(beta_2 | beta_1)\n")
cat("  H0: beta_2 = 0  (las variables omitidas no aportan)\n")
cat("  H1: beta_2 != 0 (al menos una variable omitida sí aporta)\n")
cat("  Estadístico: lambda(beta_1) - lambda(beta)\n")
cat("  Distribución bajo H0: chi^2 con r grados de libertad\n")
cat("  Criterio: RECHAZAR H0 si lambda(beta_2|beta_1) >= chi^2_{alpha, r}\n")
cat("  Modelo completo beta: base sin COT (27 variables)\n")
cat("--------------------------------------------------------------\n")

prueba2 <- function(modelo_reducido, modelo_completo, nombre, alpha = ALPHA) {
  lambda_red  <- modelo_reducido$deviance     # lambda(beta_1)
  lambda_comp <- modelo_completo$deviance     # lambda(beta)
  diff_lambda <- lambda_red - lambda_comp     # lambda(beta_2 | beta_1)
  
  p_red  <- length(coef(modelo_reducido))
  p_comp <- length(coef(modelo_completo))
  r      <- p_comp - p_red                    # parámetros omitidos
  
  critico <- qchisq(1 - alpha, df = r)
  pval    <- 1 - pchisq(diff_lambda, df = r)
  decision <- if (diff_lambda >= critico)
    "Rechazar H0 (las omitidas SI aportan)"
  else
    "No rechazar H0 (las omitidas no aportan; reducido adecuado)"
  
  cat(sprintf("\n  %s\n", nombre))
  cat(sprintf("    lambda(beta_1)            = %.3f\n", lambda_red))
  cat(sprintf("    lambda(beta)              = %.3f\n", lambda_comp))
  cat(sprintf("    lambda(beta_2 | beta_1)   = %.3f\n", diff_lambda))
  cat(sprintf("    r (parametros omitidos)   = %d - %d = %d\n", p_comp, p_red, r))
  cat(sprintf("    chi^2 critico (a=%.2f)    = %.3f\n", alpha, critico))
  cat(sprintf("    p-valor                   = %s\n",
              format.pval(pval, digits = 4, eps = 1e-4)))
  cat(sprintf("    Decision                  = %s\n", decision))
  
  data.frame(Modelo = nombre,
             diff_lambda = round(diff_lambda, 3),
             r = r, critico = round(critico, 3),
             p_valor = pval, Decision = decision,
             stringsAsFactors = FALSE)
}

p2_m1  <- prueba2(modelo_1,        modelo_base_sinCOT, "Modelo 1: Teórico")
p2_m2s <- prueba2(modelo_2_sinCOT, modelo_base_sinCOT, "Modelo 2 sin COT (AIC)")
p2_m3s <- prueba2(modelo_3_sinCOT, modelo_base_sinCOT, "Modelo 3 sin COT (BIC)")

# ------------------------------------------------------------------------------
# TABLAS RESUMEN
# ------------------------------------------------------------------------------
cat("\n\n--------------------------------------------------------------\n")
cat("RESUMEN PRUEBA 1 (desviación del modelo)\n")
cat("--------------------------------------------------------------\n")
res_p1 <- rbind(p1_m1, p1_m2s, p1_m3s)
res_p1$p_valor <- format.pval(res_p1$p_valor, digits = 4, eps = 1e-4)
print(res_p1, row.names = FALSE)

cat("\n--------------------------------------------------------------\n")
cat("RESUMEN PRUEBA 2 (subconjuntos vs base sin COT)\n")
cat("--------------------------------------------------------------\n")
res_p2 <- rbind(p2_m1, p2_m2s, p2_m3s)
res_p2$p_valor <- format.pval(res_p2$p_valor, digits = 4, eps = 1e-4)
print(res_p2, row.names = FALSE)

# ==============================================================================
# 9: PSEUDO-R^2 (McFadden, Cox-Snell, Nagelkerke)
# ==============================================================================
# En regresión logística, NO existe un único R^2 con todas las propiedades
# del coeficiente de determinación de regresión lineal. Las notas de clase
# describen 3 alternativas, las tres se calculan a continuación.
#
# 
#
#   McFadden (1974, 1977):
#       R2_McF = 1 - ln(L_modelo) / ln(L_nulo)
#       Equivalente: R2_McF = 1 - deviance_modelo / deviance_nulo
#       Rango: [0, 1) - nunca alcanza 1 en datos reales.
#
#   Cox & Snell (1989):
#       R2_CS = 1 - (L_nulo / L_modelo)^(2/n)
#       Techo teórico: R2_CS_max = 1 - L_nulo^(2/n)
#       Para datos binarios desbalanceados, el máximo está bien debajo de 1.
#
#   Nagelkerke (1991):
#       R2_N = R2_CS / R2_CS_max
#       Es Cox-Snell reescalado para tener rango [0, 1].
#
# Interpretación:
#   - McFadden: McFadden (1977, Cowles Discussion Paper 474) escribió
#     "values from .2 to .4 represent excellent fit". Esta es la única
#     guía publicada por un autor original.
#   - Cox-Snell y Nagelkerke: Hosmer, Lemeshow & Sturdivant (2013, cap 5)
#     advierten contra interpretar valores absolutos. Su uso legítimo es
#     COMPARAR modelos sobre el mismo conjunto de datos.
# ==============================================================================

cat("\n==============================================================\n")
cat("9: PSEUDO-R^2 (McFadden, Cox-Snell, Nagelkerke)\n")
cat("==============================================================\n")

# Modelo nulo (solo intercepto): denominador común para todos los pseudo-R^2
modelo_nulo <- glm(y ~ 1, data = train, family = binomial(link = "logit"))
logL_nulo   <- as.numeric(logLik(modelo_nulo))
deviance_nulo <- modelo_nulo$deviance

cat(sprintf("\n  Modelo nulo (solo intercepto):\n"))
cat(sprintf("    log-verosimilitud  = %.3f\n", logL_nulo))
cat(sprintf("    deviance           = %.3f\n", deviance_nulo))

# Función que calcula los 3 pseudo-R^2 para un modelo dado
pseudo_r2 <- function(modelo, nombre, n) {
  
  logL_mod  <- as.numeric(logLik(modelo))
  dev_mod   <- modelo$deviance
  
  # McFadden
  r2_mcf <- 1 - dev_mod / deviance_nulo
  
  # Cox-Snell
  r2_cs  <- 1 - exp((2/n) * (logL_nulo - logL_mod))
  
  # Cox-Snell máximo teórico (techo de R2_CS para esta muestra)
  r2_cs_max <- 1 - exp((2/n) * logL_nulo)
  
  # Nagelkerke
  r2_n <- r2_cs / r2_cs_max
  
  # Interpretación de McFadden según McFadden (1977):
  #   < 0.2  -> ajuste regular
  #   0.2-0.4 -> excelente
  #   > 0.4  -> excelente y poco común
  interp_mcf <- if (r2_mcf < 0.2) "regular"
  else if (r2_mcf <= 0.4) "excelente (McFadden 1977)"
  else "muy alto"
  
  cat(sprintf("\n  %s\n", nombre))
  cat(sprintf("    log-verosimilitud   = %.3f\n", logL_mod))
  cat(sprintf("    deviance            = %.3f\n", dev_mod))
  cat(sprintf("    R2 McFadden         = %.4f   [%s]\n", r2_mcf, interp_mcf))
  cat(sprintf("    R2 Cox-Snell        = %.4f   (max teorico: %.4f)\n",
              r2_cs, r2_cs_max))
  cat(sprintf("    R2 Nagelkerke       = %.4f\n", r2_n))
  
  data.frame(
    Modelo      = nombre,
    R2_McFadden = round(r2_mcf, 4),
    R2_CoxSnell = round(r2_cs, 4),
    R2_Nagelkerke = round(r2_n, 4),
    stringsAsFactors = FALSE
  )
}

# Aplicamos a los 3 modelos finales
pr2_m1  <- pseudo_r2(modelo_1,        "Modelo 1: Teórico",      n_train)
pr2_m2s <- pseudo_r2(modelo_2_sinCOT, "Modelo 2 sin COT (AIC)", n_train)
pr2_m3s <- pseudo_r2(modelo_3_sinCOT, "Modelo 3 sin COT (BIC)", n_train)

# ------------------------------------------------------------------------------
# Tabla comparativa final
# ------------------------------------------------------------------------------
cat("\n--------------------------------------------------------------\n")
cat("RESUMEN: PSEUDO-R^2 DE LOS 3 MODELOS FINALES\n")
cat("--------------------------------------------------------------\n")
cat("Nota: McFadden (1977) sugiere 0.2-0.4 = ajuste excelente.\n")
cat("Cox-Snell y Nagelkerke deben usarse para COMPARAR modelos sobre\n")
cat("el mismo dataset, no para interpretar valores absolutos\n")
cat("(Hosmer, Lemeshow & Sturdivant 2013).\n\n")

resumen_r2 <- rbind(pr2_m1, pr2_m2s, pr2_m3s)
print(resumen_r2, row.names = FALSE)

# Identificar el "ganador" en cada métrica
cat("\n  Mayor R^2 McFadden:    ",
    resumen_r2$Modelo[which.max(resumen_r2$R2_McFadden)], "\n")
cat("  Mayor R^2 Cox-Snell:   ",
    resumen_r2$Modelo[which.max(resumen_r2$R2_CoxSnell)], "\n")
cat("  Mayor R^2 Nagelkerke:  ",
    resumen_r2$Modelo[which.max(resumen_r2$R2_Nagelkerke)], "\n")


# ==============================================================================
# 10: VALIDACIÓN EN EL CONJUNTO DE PRUEBA (matriz de confusión)
# ==============================================================================
# El conjunto test (20%, n=484) ha permanecido congelado desde el inicio.
# Los modelos nunca lo vieron durante el entrenamiento, selección de
# variables, ni stepwise. Esto evita data leakage y nos da un estimador
# honesto del rendimiento predictivo de cada modelo en datos nuevos.
#
# Procedimiento:
#   1. Predecir P(y=1|x) en el test con cada modelo.
#   2. Convertir a clase binaria con umbral c = 0.5 (estándar).
#   3. Construir matriz de confusión 2x2.
#   4. Calcular las 10 métricas de las notas + kappa de Cohen.
#
# Las métricas siguen la nomenclatura de las notas:
#   a = TP (verdadero positivo): y=1, y_hat=1
#   b = FP (falso positivo):     y=0, y_hat=1
#   c = FN (falso negativo):     y=1, y_hat=0
#   d = TN (verdadero negativo): y=0, y_hat=0
#
# Nota sobre kappa: las notas presentan P_E con una multiplicación de los
# 4 marginales en el numerador, pero la fórmula correcta de Cohen (1960)
# es la SUMA de dos productos. Usamos la fórmula correcta.
# ==============================================================================

cat("\n==============================================================\n")
cat("10: VALIDACIÓN EN EL CONJUNTO DE PRUEBA\n")
cat("==============================================================\n")

UMBRAL_CLASIF <- 0.5
cat(sprintf("\n  Umbral de clasificación: c = %.2f (estándar)\n", UMBRAL_CLASIF))
cat(sprintf("  Tamaño del test: n = %d\n", nrow(test)))
cat(sprintf("  Prevalencia y=1 en test: %.3f\n", mean(test$y)))

# ------------------------------------------------------------------------------
# Función que evalúa un modelo en el test y devuelve todas las métricas
# ------------------------------------------------------------------------------
evaluar_modelo <- function(modelo, datos_test, nombre, umbral = UMBRAL_CLASIF) {
  
  # 1. Predicciones de probabilidad
  prob <- predict(modelo, newdata = datos_test, type = "response")
  
  # 2. Clase predicha
  y_hat <- as.integer(prob >= umbral)
  y_obs <- datos_test$y
  
  # 3. Matriz de confusión 2x2 (a, b, c, d según las notas)
  a <- sum(y_obs == 1 & y_hat == 1)   # TP
  b <- sum(y_obs == 0 & y_hat == 1)   # FP
  c <- sum(y_obs == 1 & y_hat == 0)   # FN
  d <- sum(y_obs == 0 & y_hat == 0)   # TN
  n <- a + b + c + d
  
  # 4. Métricas (10 de las notas + kappa)
  exactitud      <- (a + d) / n
  sensibilidad   <- a / (a + c)
  especificidad  <- d / (b + d)
  tasa_fp        <- b / (b + d)
  vpp            <- a / (a + b)
  vpn            <- d / (c + d)
  prevalencia    <- (a + c) / n
  exact_balanc   <- (sensibilidad + especificidad) / 2
  tasa_deteccion <- a / n
  preval_detec   <- (a + b) / n
  
  # Kappa de Cohen (formula correcta: suma de productos en el numerador)
  pe <- ((a + c) * (a + b) + (b + d) * (c + d)) / (n^2)
  kappa <- (exactitud - pe) / (1 - pe)
  
  # Interpretación de kappa segun las notas
  interp_kappa <- if (kappa < 0)         "sin acuerdo"
  else if (kappa <= 0.2) "insignificante"
  else if (kappa <= 0.4) "discreto"
  else if (kappa <= 0.6) "moderado"
  else if (kappa <= 0.8) "sustancial"
  else                   "casi perfecto"
  
  # ----- imprimir matriz y métricas -----
  cat(sprintf("\n--- %s ---\n", nombre))
  cat("\n  Matriz de confusión (filas = observado, columnas = predicho):\n")
  m <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE,
              dimnames = list(c("y=1", "y=0"), c("y_hat=1", "y_hat=0")))
  print(m)
  cat(sprintf("\n  Total: n = %d  (a=%d, b=%d, c=%d, d=%d)\n", n, a, b, c, d))
  
  cat("\n  Métricas:\n")
  cat(sprintf("    Exactitud (accuracy)        = %.4f\n", exactitud))
  cat(sprintf("    Sensibilidad (recall)       = %.4f\n", sensibilidad))
  cat(sprintf("    Especificidad               = %.4f\n", especificidad))
  cat(sprintf("    Tasa falsos positivos       = %.4f\n", tasa_fp))
  cat(sprintf("    VPP (precision)             = %.4f\n", vpp))
  cat(sprintf("    VPN                         = %.4f\n", vpn))
  cat(sprintf("    Prevalencia (observada)     = %.4f\n", prevalencia))
  cat(sprintf("    Exactitud balanceada        = %.4f\n", exact_balanc))
  cat(sprintf("    Tasa de detección           = %.4f\n", tasa_deteccion))
  cat(sprintf("    Prevalencia de detección    = %.4f\n", preval_detec))
  cat(sprintf("    Kappa de Cohen              = %.4f   [%s]\n",
              kappa, interp_kappa))
  
  data.frame(
    Modelo                 = nombre,
    a_TP                   = a,
    b_FP                   = b,
    c_FN                   = c,
    d_TN                   = d,
    Exactitud              = round(exactitud, 4),
    Sensibilidad           = round(sensibilidad, 4),
    Especificidad          = round(especificidad, 4),
    Tasa_FP                = round(tasa_fp, 4),
    VPP                    = round(vpp, 4),
    VPN                    = round(vpn, 4),
    Exact_Balanceada       = round(exact_balanc, 4),
    Tasa_Deteccion         = round(tasa_deteccion, 4),
    Prev_Deteccion         = round(preval_detec, 4),
    Kappa                  = round(kappa, 4),
    Kappa_Interp           = interp_kappa,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Aplicar a los 3 modelos finales
# ------------------------------------------------------------------------------
ev_m1  <- evaluar_modelo(modelo_1,        test, "Modelo 1: Teórico")
ev_m2s <- evaluar_modelo(modelo_2_sinCOT, test, "Modelo 2 sin COT (AIC)")
ev_m3s <- evaluar_modelo(modelo_3_sinCOT, test, "Modelo 3 sin COT (BIC)")

# ------------------------------------------------------------------------------
# Tabla resumen comparativa
# ------------------------------------------------------------------------------
cat("\n--------------------------------------------------------------\n")
cat("RESUMEN COMPARATIVO EN TEST\n")
cat("--------------------------------------------------------------\n")
resumen_test <- rbind(ev_m1, ev_m2s, ev_m3s)

# columnas clave para comparar
print(resumen_test[, c("Modelo", "Exactitud", "Sensibilidad", "Especificidad",
                       "Exact_Balanceada", "Kappa", "Kappa_Interp")],
      row.names = FALSE)

# Identificar ganadores
cat("\n  Mayor exactitud:           ",
    resumen_test$Modelo[which.max(resumen_test$Exactitud)], "\n")
cat("  Mayor sensibilidad:        ",
    resumen_test$Modelo[which.max(resumen_test$Sensibilidad)], "\n")
cat("  Mayor especificidad:       ",
    resumen_test$Modelo[which.max(resumen_test$Especificidad)], "\n")
cat("  Mayor exactitud balanceada:",
    resumen_test$Modelo[which.max(resumen_test$Exact_Balanceada)], "\n")
cat("  Mayor kappa:               ",
    resumen_test$Modelo[which.max(resumen_test$Kappa)], "\n")

# ==============================================================================
# 11: VERIFICACIÓN DE SUPUESTOS DE LOS MODELOS FINALES
# ==============================================================================
# La regresión logística tiene MENOS supuestos que la regresión lineal
# clásica porque no asume normalidad de errores ni homocedasticidad.
# Los supuestos a verificar son:
#
#   (A) LINEALIDAD DEL LOGIT: la relación entre cada variable predictora
#       CONTINUA y el log-odds de y debe ser lineal. Se verifica con la
#       prueba de Box-Tidwell, que añade el término x*ln(x) al modelo y
#       prueba si su coeficiente es significativo. Si lo es, hay no
#       linealidad. Se complementa con gráficas de logit empírico.
#
#   (B) AUSENCIA DE MULTICOLINEALIDAD: las regresoras no deben estar
#       fuertemente correlacionadas entre sí. Se verifica con VIF.
#
#   (C) INDEPENDENCIA: las observaciones deben ser independientes entre
#       sí. Se verifica con un análisis exploratorio de mediciones
#       repetidas por sitio.
#
# Las binarias NO se evalúan en linealidad del logit (no aplica con 2
# niveles).
# ==============================================================================

cat("\n==============================================================\n")
cat("11: VERIFICACIÓN DE SUPUESTOS\n")
cat("==============================================================\n")

library(car)  # boxTidwell, vif

# ------------------------------------------------------------------------------
# 11A: LINEALIDAD DEL LOGIT (Box-Tidwell)
# ------------------------------------------------------------------------------
# Box-Tidwell prueba: H0: la variable cumple linealidad del logit
# Estadístico: significancia del término x*ln(x) en el modelo extendido.
#   - p < 0.05 -> rechazar H0 -> hay no linealidad
#   - p >= 0.05 -> no rechazar H0 -> linealidad se cumple
#
# Implementación manual: para cada variable continua x_j de un modelo,
# ajustamos el modelo extendido: y ~ ... + x_j + x_j:log(x_j).
# Si el coeficiente de x_j:log(x_j) es significativo (Wald), hay no
# linealidad en x_j.
#
# Detalle técnico: log(0) no existe. Para variables con ceros sumamos un
# epsilon pequeño antes de tomar log. Esto introduce un sesgo mínimo.

cat("\n--------------------------------------------------------------\n")
cat("11A: LINEALIDAD DEL LOGIT (Box-Tidwell)\n")
cat("--------------------------------------------------------------\n")
cat("  H0: linealidad del logit  vs  H1: no linealidad\n")
cat("  Decision: si p_valor del termino x*log(x) < 0.05 -> rechazar H0\n\n")

box_tidwell_modelo <- function(modelo, datos, nombre) {
  
  cat(sprintf("\n--- %s ---\n", nombre))
  
  # variables del modelo
  vars_modelo <- attr(terms(modelo), "term.labels")
  vars_modelo <- gsub("`", "", vars_modelo)
  
  # filtrar solo las continuas (las binarias 0/1 no aplican)
  vars_continuas <- vars_modelo[
    sapply(vars_modelo, function(v) !all(datos[[v]] %in% c(0, 1)))
  ]
  
  resultados <- data.frame(
    Variable = character(), p_valor = numeric(),
    Decision = character(), stringsAsFactors = FALSE
  )
  
  formula_orig <- formula(modelo)
  
  for (v in vars_continuas) {
    # epsilon para variables con ceros
    valores <- datos[[v]]
    epsilon <- if (min(valores) <= 0) 1e-6 else 0
    
    # crear término de interacción x*log(x) en una columna nueva
    nombre_int <- paste0(v, "_xlogx")
    datos_aux  <- datos
    datos_aux[[nombre_int]] <- valores * log(valores + epsilon)
    
    # modelo extendido: y ~ todas las del modelo + x*log(x)
    formula_ext <- update(formula_orig,
                          as.formula(paste(". ~ . +", paste0("`", nombre_int, "`"))))
    
    modelo_ext <- tryCatch(
      glm(formula_ext, data = datos_aux,
          family = binomial(link = "logit")),
      error = function(e) NULL,
      warning = function(w) {
        glm(formula_ext, data = datos_aux,
            family = binomial(link = "logit"))
      }
    )
    
    if (is.null(modelo_ext)) {
      next
    }
    
    # extraer p-valor del termino x*log(x)
    coefs_ext <- summary(modelo_ext)$coefficients
    fila_int  <- coefs_ext[grepl(nombre_int, rownames(coefs_ext), fixed = TRUE), ,
                           drop = FALSE]
    
    if (nrow(fila_int) == 0) next
    
    pv <- fila_int[1, "Pr(>|z|)"]
    decision <- if (pv < 0.05) "Rechazar H0 (no linealidad)" else "No rechazar H0 (lineal OK)"
    
    resultados <- rbind(resultados, data.frame(
      Variable = v,
      p_valor  = pv,
      Decision = decision,
      stringsAsFactors = FALSE
    ))
  }
  
  resultados$p_valor <- format.pval(resultados$p_valor, digits = 4, eps = 1e-4)
  print(resultados, row.names = FALSE)
  
  invisible(resultados)
}

bt_m1  <- box_tidwell_modelo(modelo_1,        train, "Modelo 1: Teórico")
bt_m2s <- box_tidwell_modelo(modelo_2_sinCOT, train, "Modelo 2 sin COT (AIC)")
bt_m3s <- box_tidwell_modelo(modelo_3_sinCOT, train, "Modelo 3 sin COT (BIC)")

# ------------------------------------------------------------------------------
# 11B: MULTICOLINEALIDAD (VIF) DENTRO DE CADA MODELO FINAL
# ------------------------------------------------------------------------------
# El VIF del pool inicial (sección 5) ya filtró redundancias globales.
# Aquí calculamos el VIF DENTRO de cada modelo final como diagnóstico
# adicional: con menos variables el VIF puede bajar todavia más.
# Criterio: VIF <= 10 (estándar de la rúbrica).

cat("\n\n--------------------------------------------------------------\n")
cat("11B: VIF DENTRO DE CADA MODELO FINAL\n")
cat("--------------------------------------------------------------\n")
cat("  Criterio: VIF <= 10 indica ausencia de multicolinealidad crítica\n\n")

vif_modelo <- function(modelo, nombre) {
  cat(sprintf("\n--- %s ---\n", nombre))
  vifs <- vif(modelo)
  vifs_ord <- sort(vifs, decreasing = TRUE)
  for (v in names(vifs_ord)) {
    flag <- if (vifs_ord[v] > 10) " ⚠ alto" else ""
    cat(sprintf("    %-25s VIF = %6.3f%s\n", gsub("`", "", v), vifs_ord[v], flag))
  }
  cat(sprintf("    -> max VIF = %.2f  |  min VIF = %.2f\n",
              max(vifs), min(vifs)))
}

vif_modelo(modelo_1,        "Modelo 1: Teórico")
vif_modelo(modelo_2_sinCOT, "Modelo 2 sin COT (AIC)")
vif_modelo(modelo_3_sinCOT, "Modelo 3 sin COT (BIC)")

# ------------------------------------------------------------------------------
# 11C: INDEPENDENCIA - análisis de mediciones repetidas por sitio
# ------------------------------------------------------------------------------
# La regresión logística asume que las observaciones son independientes
# entre sí. Si la base contiene mediciones repetidas del mismo sitio
# (en distintos años), esas observaciones NO son estrictamente
# independientes: tienden a parecerse por características del sitio.
#
# Si la mayoría de sitios tiene UNA medición -> independencia aprox OK.
# Si hay muchos sitios con varias mediciones -> violación; sería
# necesario un modelo de efectos mixtos (no cubierto en este curso).
# Se reporta como limitación del modelo.

cat("\n\n--------------------------------------------------------------\n")
cat("11C: INDEPENDENCIA (mediciones repetidas por sitio)\n")
cat("--------------------------------------------------------------\n")

# contar mediciones por sitio en el train
conteo_sitios <- table(train[["CLAVE SITIO"]])

cat(sprintf("  Total observaciones en train: %d\n", nrow(train)))
cat(sprintf("  Sitios únicos: %d\n", length(conteo_sitios)))
cat(sprintf("  Mediciones por sitio:\n"))
cat(sprintf("    media:   %.2f\n", mean(conteo_sitios)))
cat(sprintf("    mediana: %.0f\n", median(conteo_sitios)))
cat(sprintf("    mín:     %d\n", min(conteo_sitios)))
cat(sprintf("    máx:     %d\n", max(conteo_sitios)))

# distribución
cat("\n  Distribución del nro de mediciones por sitio:\n")
dist_med <- table(conteo_sitios)
for (k in names(dist_med)) {
  cat(sprintf("    %2s mediciones: %4d sitios (%.1f%% de los sitios)\n",
              k, dist_med[k], 100 * dist_med[k] / length(conteo_sitios)))
}

# proporción de observaciones que vienen de sitios con N>1
sitios_repetidos <- names(conteo_sitios[conteo_sitios > 1])
n_obs_repetidas  <- sum(train[["CLAVE SITIO"]] %in% sitios_repetidos)
prop_repetidas   <- n_obs_repetidas / nrow(train)

cat(sprintf("\n  Observaciones de sitios con >1 medición: %d / %d (%.1f%%)\n",
            n_obs_repetidas, nrow(train), 100 * prop_repetidas))

# Diagnóstico
cat("\n  Diagnóstico:\n")
if (prop_repetidas < 0.20) {
  cat("    -> Independencia se cumple aproximadamente.\n")
} else if (prop_repetidas < 0.50) {
  cat("    -> Hay dependencia moderada por sitio. Reportar como limitación.\n")
} else {
  cat("    -> Dependencia fuerte. Idealmente modelo de efectos mixtos.\n")
  cat("    -> Se reporta como limitación del estudio.\n")
}