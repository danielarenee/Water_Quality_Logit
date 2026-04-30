# ==============================================================================
# PROYECTO DE REGRESIÓN LOGÍSTICA
# DANIELA RENÉE Y MARIO SÁNCHEZ
# ECONOMETRÍA 1
# ==============================================================================
#Ruta Dani Abajo
#setwd("/Users/danielarenee/Desktop/Water_Quality_Logit")
#Ruta Mario Abajo
setwd("C:/Users/msfcl/OneDrive/Escritorio/Water_Quality_Logit")

library(car)      
library(corrplot)
library(dplyr)
library(pscl)

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


# 6. AJUSTE DE LOS 4 MODELOS CANDIDATOS

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
cat(sprintf("  Variables: %d | AIC: %.2f | BIC: %.2f\n", length(vars_m1), AIC(modelo_1), BIC(modelo_1)))

# VIF para multicolinealidad 
print(car::vif(modelo_1)) 

# ------------------------------------------------------------------------------
# MODELO 2: Stepwise por AIC (bidireccional)
# ------------------------------------------------------------------------------
# step() arranca desde modelo_base y en cada paso considera agregar o quitar
cat("\n--- Modelo 2: Stepwise AIC (bidireccional) ---\n")

modelo_2 <- step(modelo_base, direction = "both", trace = 0, k = 2)
vars_m2  <- gsub("`", "", names(coef(modelo_2))[-1]) 
cat(sprintf("  Variables: %d | AIC: %.2f | BIC: %.2f\n", length(vars_m2), AIC(modelo_2), BIC(modelo_2)))

# VIF para multicolinealidad
print(car::vif(modelo_2))

# Chi cuadrado para independencia de categoricas
vars_cat2 <- c("COT_NOM", "TOX_NOM", "TIPO_LOTICO")
combi <- combn(vars_cat2, 2)

for(i in 1:ncol(combi)) {
  v1 <- combi[1,i]; v2 <- combi[2,i]
  cat(sprintf("\n> %s vs %s:\n", v1, v2))
  print(chisq.test(table(train[[v1]], train[[v2]])))
}

# Hay redundancia CRITICA entre COT_NOM y TOX_NOM
# independencia entre tipo_lotico y las otras 2 

# ------------------------------------------------------------------------------
# MODELO 3: Stepwise por BIC (bidireccional)
# ------------------------------------------------------------------------------
# BIC penaliza más la complejidad que AIC
cat("\n--- Modelo 3: Stepwise BIC (bidireccional) ---\n")

n_train  <- nrow(train)
modelo_3 <- step(modelo_base, direction = "both", trace = 0, k = log(n_train))
vars_m3  <- gsub("`", "", names(coef(modelo_3))[-1])

cat(sprintf("  Variables: %d | AIC: %.2f | BIC: %.2f\n", length(vars_m3), AIC(modelo_3), BIC(modelo_3)))
cat("  Variables:", paste(vars_m3, collapse = ", "), "\n")

# VIF para multicolinealidad
print(car::vif(modelo_3))

# Chi cuadrado para independencia (COT_NOM vs TIPO_LOTICO)
print(chisq.test(table(train[["COT_NOM"]], train[["TIPO_LOTICO"]])))
# independencia


# 7. INFERENCIA DE WALD - PRUEBA DE SIGNIFICANCIA INDIVIDUAL

# Para cada coeficiente del modelo se prueba:
#   H0: beta_j = 0    vs    H1: beta_j != 0
#   Z_0 = beta_hat_j / se(beta_hat_j)
#   alpha 0.05

cat("7: INFERENCIA DE WALD POR COEFICIENTE\n")

# función para construir la tabla de Wald de cualquier glm
tabla_wald <- function(modelo, nombre_modelo, alpha = 0.05) {
  
  coefs <- summary(modelo)$coefficients
  tabla <- data.frame(
    Variable     = rownames(coefs),
    Beta         = round(coefs[, "Estimate"],   4),
    SE           = round(coefs[, "Std. Error"], 4),
    Z            = round(coefs[, "z value"],    3),
    p_valor      = coefs[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
    tabla$Decision <- ifelse(tabla$p_valor < alpha,
                           "Rechazar H0", "No rechazar H0")
  
  tabla$Sig <- cut(tabla$p_valor,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns"),
                   right  = FALSE)
  
    tabla$p_valor <- format.pval(tabla$p_valor, digits = 4, eps = 1e-4)
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

# aquí notamos que COT_NOM es cuasi-completa, y está incluida en M2 y M3
# esto porque el Carbono Orgánico Total (COT) y la Demanda Química de Oxígeno 
# son medidas del mismo fenómeno: la carga orgánica del agua. 

# Construimos versiones "sin COT" de M2 y M3 p

cat(": VERSIONES SIN COT_NOM\n")

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

# comparamos AIC para cuantificar el impacto de cot_nom y justificar 
# su eliminación de aquí en adelante

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

cat("\n  Penalización por quitar COT_NOM (en AIC):\n")
cat(sprintf("    M2 -> M2 sin COT:  delta AIC = %+.2f\n",
            AIC(modelo_2_sinCOT) - AIC(modelo_2)))
cat(sprintf("    M3 -> M3 sin COT:  delta AIC = %+.2f\n",
            AIC(modelo_3_sinCOT) - AIC(modelo_3)))

# 8: PRUEBAS DE DESVIACIÓN 

# PRUEBA 1: Desviación del modelo (vs modelo saturado teórico):
#   Bajo H0: el modelo ajustado es adecuado
#   Criterio: NO RECHAZAR si lambda(beta) <= chi^2_{alpha, n-p}

# PRUEBA 2 - Subconjuntos de parámetros:
#   H0: beta_2 = 0  vs  H1: beta_2 != 0
#   lambda(beta_2 | beta_1) = lambda(beta_1) - lambda(beta)
#   Bajo H0, distribuye chi^2_r donde r = nro de parámetros omitidos
#   Criterio: RECHAZAR H0 si lambda(beta_2 | beta_1) >= chi^2_{alpha, r}

# El "modelo completo" beta para la Prueba 2 es el modelo de 27 variables

cat("8: PRUEBAS DE DESVIACIÓN\n")
ALPHA <- 0.05

# Modelo base sin COT_NOM (27 variables)

modelo_base_sinCOT <- glm(y ~ . - COT_NOM, 
                          data = train[, c("y", vars_actuales)], 
                          family = binomial)
cat(sprintf("\nDeviance: %.2f   AIC: %.2f\n", deviance(modelo_base_sinCOT), AIC(modelo_base_sinCOT)))

# ------------------------------------------------------------------------------
# PRUEBA 1: Desviación del modelo (vs saturado teórico)
# ------------------------------------------------------------------------------

prueba1 <- function(modelo, nombre, n, alpha = ALPHA) {
  lambda  <- modelo$deviance
  p       <- length(coef(modelo))      # parámetros (incluye intercepto)
  gl      <- n - p
  critico <- qchisq(1 - alpha, df = gl)
  pval    <- 1 - pchisq(lambda, df = gl)
  decision <- if (lambda <= critico) "No rechazar H0 (modelo adecuado)"
  else                   "Rechazar H0 (modelo NO adecuado)"
  
  cat(sprintf("\n  %s\n", nombre))
  cat(sprintf("    Decision            = %s\n", decision))
}

n_train <- nrow(train)
p1_m1  <- prueba1(modelo_1,        "Modelo 1: Teórico",        n_train)
p1_m2s <- prueba1(modelo_2_sinCOT, "Modelo 2 (AIC)",   n_train)
p1_m3s <- prueba1(modelo_3_sinCOT, "Modelo 3 (BIC)",   n_train)

# ------------------------------------------------------------------------------
# PRUEBA 2: Subconjuntos de parámetros (vs modelo base sin COT)
# ------------------------------------------------------------------------------

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
  cat(sprintf("    Decision                  = %s\n", decision))
}

p2_m1  <- prueba2(modelo_1,        modelo_base_sinCOT, "Modelo 1: Teórico")
p2_m2s <- prueba2(modelo_2_sinCOT, modelo_base_sinCOT, "Modelo 2 sin COT (AIC)")
p2_m3s <- prueba2(modelo_3_sinCOT, modelo_base_sinCOT, "Modelo 3 sin COT (BIC)")

# 9: PSEUDO-R^2 

# McFadden: R2_McF = 1 - ln(L_modelo) / ln(L_nulo)
# Cox & Snell (r2ML): R2_CS = 1 - (L_nulo / L_modelo)^(2/n)
# Nagelkerke (r2CU): R2_N = R2_CS / R2_CS_max

cat("9: PSEUDO-R^2 (McFadden, Cox-Snell, Nagelkerke)\n")

pR2(modelo_1)
pR2(modelo_2_sinCOT) # ganador en las 3
pR2(modelo_3_sinCOT)

# 10: VALIDACIÓN EN EL CONJUNTO DE PRUEBA (matriz de confusión)
#
#   a = TP : y=1, y_hat=1
#   b = FP : y=0, y_hat=1 # costo economico
#   c = FN : y=1, y_hat=0 # riesgo sanitario
#   d = TN : y=0, y_hat=0

cat("10: VALIDACIÓN EN EL CONJUNTO DE PRUEBA\n")

UMBRAL_CLASIF <- 0.5

# Función que evalúa un modelo en el test y devuelve todas las métricas
evaluar_modelo <- function(modelo, datos_test, nombre, umbral = UMBRAL_CLASIF) {
  
  # 1. Predicciones de probabilidad
  prob <- predict(modelo, newdata = datos_test, type = "response")
  # 2. Clase predicha
  y_hat <- as.integer(prob >= umbral)
  y_obs <- datos_test$y
  
  # 3. Matriz de confusión 2x2 
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
  
  # Kappa de Cohen 
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
}

# Aplicar a los 3 modelos finales
ev_m1  <- evaluar_modelo(modelo_1,        test, "Modelo 1: Teórico")
ev_m2s <- evaluar_modelo(modelo_2_sinCOT, test, "Modelo 2 sin COT (AIC)")
ev_m3s <- evaluar_modelo(modelo_3_sinCOT, test, "Modelo 3 sin COT (BIC)")


# 11: VERIFICACIÓN DE SUPUESTOS DE LOS MODELOS FINALES

cat("11: VERIFICACIÓN DE SUPUESTOS\n")

# MULTICOLINEALIDAD
print(car::vif(modelo_1))
print(car::vif(modelo_2_sinCOT))
print(car::vif(modelo_3_sinCOT))

# LINEALIDAD (Box-Tidwell)

box_tidwell_modelo <- function(modelo, datos, nombre) {
  cat(sprintf("\n--- %s ---\n", nombre))
  
  # 1. Definir variables continuas (excluyendo las 3 categóricas y la respuesta)
  vars_modelo <- gsub("`", "", attr(terms(modelo), "term.labels"))
  mis_categoricas <- c("COT_NOM", "TOX_NOM", "TIPO_LOTICO")
  vars_cont <- setdiff(vars_modelo, mis_categoricas)
  
  # 2. Crear términos x*log(x+1) dinámicamente
  # Usamos log(x+1) para evitar problemas con ceros 
  terminos_bt <- sprintf("I(`%s` * log(`%s` + 1))", vars_cont, vars_cont)
  
  # 3. Construir y ajustar el modelo extendido
  formula_ext <- update(formula(modelo), paste(". ~ . +", paste(terminos_bt, collapse = " + ")))
  modelo_ext <- glm(formula_ext, data = datos, family = binomial)
  coefs <- summary(modelo_ext)$coefficients
  
  # 4. Filtrar y presentar resultados
  # Buscamos las filas que contienen "log" en el nombre del coeficiente
  filas_bt <- coefs[grepl("log", rownames(coefs)), , drop = FALSE]
  
  resultados <- data.frame(
    Variable = vars_cont,
    p_valor  = round(filas_bt[, "Pr(>|z|)"], 4),
    Decision = ifelse(filas_bt[, "Pr(>|z|)"] < 0.05, "Rechazar H0 (No lineal)", "Lineal OK"),
    stringsAsFactors = FALSE
  )
  
  print(resultados, row.names = FALSE)
  return(invisible(resultados))
}

# Ejecución
bt_m1  <- box_tidwell_modelo(modelo_1, train, "Modelo 1: Teórico")
bt_m2s <- box_tidwell_modelo(modelo_2_sinCOT, train, "Modelo 2: AIC")
bt_m3s <- box_tidwell_modelo(modelo_3_sinCOT, train, "Modelo 3: BIC")

# Notamos que no todos cumplen el supuesto de linealidad...

# INVESTIGACIÓN DE NO-LINEALIDAD CON LOG-TRANSFORMACIÓN
# Reespecificamos el modelo aplicando log(x+1) a las variables que violan el supuesto
# elegimos esta porque es interpretable y no agrega parámetros al modelo 

cat("11A: INVESTIGACIÓN DE NO-LINEALIDAD CON LOG(X+1)\n")

vars_a_transformar <- c("SST", "NI_TOT", "P_TOT", "ORTO_PO4", "ABS_UV", "CONDUC_CAMPO", "COLOR_VER", "OD_mgL")

for(v in vars_a_transformar) {
  if(v %in% names(train)) {
    train[[paste0("LOG_", v)]] <- log1p(train[[v]])
    test[[paste0("LOG_", v)]]  <- log1p(test[[v]])
  }
}

# Qué onda con SST?...
summary(train$LOG_SST)
outliers_log <- train$LOG_SST[train$LOG_SST > mean(train$LOG_SST) + 3*sd(train$LOG_SST)]
outliers_log
boxplot(train$LOG_SST, col = "tomato", main = "Detección de Outliers en SST")

# hay 2 outliers que superan las 3 desviaciones estandar
# los identificamos y eliminamos
train_limpio <- train[train$LOG_SST < 9, ]
cat(sprintf("Observaciones eliminadas: %d\n", nrow(train) - nrow(train_limpio)))

# M1: Teórico
m1_final <- glm(y ~ TEMP_AGUA + LOG_SST + LOG_COLOR_VER + pH_CAMPO + 
                  LOG_CONDUC_CAMPO + LOG_OD_mgL + NI_TOT + E_COLI, 
                data = train_limpio, family = binomial)

# M2: Stepwise AIC 
m2_final <- glm(y ~ N_NH3 + N_NO2 + LOG_P_TOT + LOG_ORTO_PO4 + LOG_ABS_UV + 
                  LOG_CONDUC_CAMPO + SAAM + LOG_SST + CD_TOT + LOG_NI_TOT + 
                  TOX_NOM + TIPO_LOTICO, 
                data = train_limpio, family = binomial)

# M3: Stepwise BIC
m3_final <- glm(y ~ N_NH3 + N_NO2 + LOG_P_TOT + LOG_ORTO_PO4 + LOG_ABS_UV + 
                  LOG_CONDUC_CAMPO + SAAM + LOG_SST + 
                  TOX_NOM + TIPO_LOTICO, 
                data = train_limpio, family = binomial)

# COMPARACIÓN DE DESEMPEÑO

comparar_modelos_fijo <- function(m_lin, m_log, nombre) {
  cat(sprintf("\n--- %s: Lineal vs Log ---\n", nombre))
  df <- data.frame(
    Metrica = c("AIC", "BIC", "Deviance"),
    Lineal  = c(AIC(m_lin), BIC(m_lin), m_lin$deviance),
    Log     = c(AIC(m_log), BIC(m_log), m_log$deviance)
  )
  df$Lineal <- round(df$Lineal, 2)
  df$Log    <- round(df$Log, 2)
  
  print(df, row.names = FALSE)
  cat(sprintf("Delta AIC: %.2f\n", AIC(m_log) - AIC(m_lin)))
}

comparar_modelos_fijo(modelo_1, m1_final, "MODELO 1")
comparar_modelos_fijo(modelo_2_sinCOT, m2_final, "MODELO 2")
comparar_modelos_fijo(modelo_3_sinCOT, m3_final, "MODELO 3")

# VERIFICACIÓN DE LINEALIDAD (Box-Tidwell Post-Log)
bt_m1_log <- box_tidwell_modelo(m1_final, train, "Modelo 1 (Log)")
bt_m2_log <- box_tidwell_modelo(m2_final, train, "Modelo 2 (Log)")
bt_m3_log <- box_tidwell_modelo(m3_final, train, "Modelo 3 (Log)")

# EVALUAR PREDICCIÓN

evaluar_prediccion <- function(modelo, nombre_modelo) {
  prob <- predict(modelo, newdata = test, type = "response")
  pred <- as.integer(prob >= 0.5)
  
  # metricas
  acc  <- mean(pred == test$y)
  sens <- sum(test$y == 1 & pred == 1) / sum(test$y == 1) # Recall
  spec <- sum(test$y == 0 & pred == 0) / sum(test$y == 0) # Especificidad
  
  cat(sprintf("\n--- Métricas en Test: %s ---\n", nombre_modelo))
  cat(sprintf("  Accuracy:      %.4f\n", acc))
  cat(sprintf("  Sensibilidad:  %.4f (Recall)\n", sens))
  cat(sprintf("  Especificidad: %.4f\n", spec))
}

evaluar_prediccion(m1_final, "Modelo 1 (Log)")
evaluar_prediccion(m2_final, "Modelo 2 (Log)")
evaluar_prediccion(m3_final, "Modelo 3 (Log)")

#checar vif
vif(m1_final)
vif(m2_final)
vif(m3_final)

# EXTRA: INDEPENDENCIA ENTRE LA VARIABLE RESPUESTA Y EL TIPO DE CUERPO DE AGUA
# Verificamos si la clasificación de calidad (DQO, base de y) depende del
# tipo de cuerpo de agua. Aplicamos sobre el dataset completo train+test

cat("EXTRA: INDEPENDENCIA ENTRE y Y TIPO_LOTICO (base completa)\n")
cat("  H0: y (semáforo DQO) y TIPO_LOTICO son independientes\n")

base_completa <- rbind(train, test)
tabla <- table(base_completa$TIPO_LOTICO, base_completa$y)
chi_test <- chisq.test(tabla, correct = FALSE)

print(round(chi_test$expected, 2)) # vemos que son >=5, chi cuadrada es valida

cat(sprintf("Estadístico: %.3f\n", chi_test$statistic))
cat(sprintf("p-valor:     %s\n", format.pval(chi_test$p.value, eps = 1e-4)))

# no se rechaza al 0.001

# ==============================================================================
# CÓDIGO OPTIMIZADO: WALD, PSEUDO R^2 Y DESVIACIÓN
# ==============================================================================

library(pscl)

# Agrupamos los modelos en una lista para iterar sin repetir código
modelos <- list(M1_Teorico = m1_final, M2_AIC = m2_final, M3_BIC = m3_final)

# ------------------------------------------------------------------------------
# TAREA A: INFERENCIA DE WALD
# summary() de R base ya calcula la matriz Hessiana y hace la prueba de Wald.
# Extraemos la tabla de coeficientes directamente.
# ------------------------------------------------------------------------------
cat("\n--- INFERENCIA DE WALD ---\n")
lapply(modelos, function(m) round(summary(m)$coefficients, 4))


# ------------------------------------------------------------------------------
# TAREA B: PSEUDO R^2
# pR2() calcula McFadden, r2ML (Cox & Snell) y r2CU (Nagelkerke) en 1 línea.
# ------------------------------------------------------------------------------
cat("\n--- PSEUDO R^2 ---\n")
# Usamos sapply y t() para que el output sea una tabla comparativa limpia
tabla_r2 <- t(sapply(modelos, pR2))
print(round(tabla_r2[, c("llh", "McFadden", "r2ML", "r2CU")], 4))


# ------------------------------------------------------------------------------
# TAREA C: PRUEBAS DE DESVIACIÓN
# ------------------------------------------------------------------------------

cat("\n--- PRUEBA 1: DESVIACIÓN VS SATURADO TEÓRICO ---\n")
# pchisq() evalúa la deviance residual (lambda) usando los grados de libertad.
# p > 0.05 indica que NO se rechaza H0 (el modelo es adecuado)
prueba1 <- sapply(modelos, function(m) {
  p_val <- 1 - pchisq(m$deviance, m$df.residual)
  c(Deviance = m$deviance, p_valor = p_val, 
    Adecuado = ifelse(p_val > 0.05, "SI", "NO"))
})
print(t(prueba1))


cat("\n--- PRUEBA 2: SUBCONJUNTOS DE PARÁMETROS VS COMPLETO TRANSFORMADO ---\n")
# 1. Definimos el modelo completo estrictamente anidado
vars_completas_log <- c("COLI_FEC", "COLI_TOT", "E_COLI", "N_NH3", "N_NO2", "N_TOTK", 
                        "TOX_D_48_UT", "LOG_P_TOT", "LOG_ORTO_PO4", "LOG_COLOR_VER", "LOG_ABS_UV", 
                        "LOG_CONDUC_CAMPO", "pH_CAMPO", "LOG_OD_mgL", "SAAM", "LOG_SST", "CD_TOT", 
                        "CR_TOT", "HG_TOT", "LOG_NI_TOT", "CN_TOT", "DUR_TOT", "TEMP_AMB", 
                        "TEMP_AGUA", "CAUDAL", "TOX_NOM", "TIPO_LOTICO")

form_completa_log <- as.formula(paste("y ~", paste(vars_completas_log, collapse = " + ")))
modelo_completo_log <- glm(form_completa_log, data = train_limpio, family = binomial)

# 2. Comparamos devianzas con anova()
# p < 0.05 indica que SE RECHAZA H0 (las variables omitidas sí aportaban)
lapply(modelos, function(m) anova(m, modelo_completo_log, test = "Chisq"))

# ==============================================================================
# EVALUACIÓN EXHAUSTIVA DE LA MATRIZ DE CONFUSIÓN (TEST)
# Cálculo de las 10 métricas de clase + Índice Kappa de Cohen
# ==============================================================================

cat("\n--- MÉTRICAS COMPLETAS DE LA MATRIZ DE CONFUSIÓN ---\n")

umbral <- 0.5
modelos <- list(M1_Teorico = m1_final, M2_AIC = m2_final, M3_BIC = m3_final)

resultados_cm <- lapply(names(modelos), function(nombre) {
  m <- modelos[[nombre]]
  
  # Predicciones
  prob_pred <- predict(m, newdata = test, type = "response")
  y_pred <- ifelse(prob_pred >= umbral, 1, 0)
  y_obs <- test$y
  
  # 1. Elementos de la matriz
  a <- sum(y_obs == 1 & y_pred == 1) # Verdaderos Positivos (VP)
  b <- sum(y_obs == 0 & y_pred == 1) # Falsos Positivos (FP)
  c <- sum(y_obs == 1 & y_pred == 0) # Falsos Negativos (FN)
  d <- sum(y_obs == 0 & y_pred == 0) # Verdaderos Negativos (VN)
  n <- a + b + c + d
  
  # 2. Las 10 métricas de las notas
  exactitud         <- (a + d) / n
  sensitividad      <- a / (a + c)
  especificidad     <- d / (b + d)
  tasa_fp           <- b / (b + d)
  vpp_precision     <- a / (a + b)
  vpn               <- d / (c + d)
  prevalencia       <- (a + c) / n
  exact_balanceada  <- (sensitividad + especificidad) / 2
  tasa_deteccion    <- a / n
  prevalencia_detec <- (a + b) / n
  
  # 3. Índice Kappa de Cohen
  pe <- ((a + c) * (a + b) + (b + d) * (c + d)) / (n^2)
  kappa <- (exactitud - pe) / (1 - pe)
  
  # Estructurar resultados
  data.frame(
    Modelo           = nombre,
    Exactitud        = round(exactitud, 4),
    Sensitividad     = round(sensitividad, 4),
    Especificidad    = round(especificidad, 4),
    Tasa_FP          = round(tasa_fp, 4),
    Precisión_VPP    = round(vpp_precision, 4),
    VPN              = round(vpn, 4),
    Prevalencia      = round(prevalencia, 4),
    Exact_Balanceada = round(exact_balanceada, 4),
    Tasa_Deteccion   = round(tasa_deteccion, 4),
    Prev_Deteccion   = round(prevalencia_detec, 4),
    Kappa            = round(kappa, 4),
    stringsAsFactors = FALSE
  )
})

# Unir y transponer para facilitar la lectura de todas las métricas
tabla_completa <- do.call(rbind, resultados_cm)
print(t(tabla_completa))

# ==============================================================================
# TABLAS DE ODDS RATIOS (COCIENTE DE POSIBILIDADES)
# Basado en la relación oddsratio = exp(Beta_hat)
# ==============================================================================

# 1. Función reutilizable para extraer coeficientes, quitar intercepto y calcular OR
generar_tabla_or <- function(modelo, nombre_modelo) {
  
  # Extraer coeficientes del modelo
  coefs <- coef(modelo)
  
  # Excluir el intercepto
  coefs <- coefs[names(coefs) != "(Intercept)"]
  
  # Calcular Odds Ratios: exp(Beta_hat)
  or_vals <- exp(coefs)
  
  # Construir el data.frame con 4 decimales
  tabla <- data.frame(
    Variable = names(coefs),
    Beta_hat = round(coefs, 4),
    OR       = round(or_vals, 4),
    stringsAsFactors = FALSE
  )
  
  # Limpiar backticks de los nombres de variables
  tabla$Variable <- gsub("`", "", tabla$Variable)
  
  # Imprimir con encabezado claro
  cat(sprintf("\n--- ODDS RATIOS: %s ---\n", nombre_modelo))
  print(tabla, row.names = FALSE)
  
  return(tabla)
}

# 2. Aplicar la función a los 3 modelos finales
tabla_m1 <- generar_tabla_or(m1_final, "Modelo 1 (Teórico)")
tabla_m2 <- generar_tabla_or(m2_final, "Modelo 2 (Stepwise AIC)")
tabla_m3 <- generar_tabla_or(m3_final, "Modelo 3 (Stepwise BIC)")


# ==============================================================================
# TABLA RESUMEN COMPARATIVA (Variables en más de un modelo)
# ==============================================================================

# Extraer todas las variables únicas presentes en los 3 modelos
todas_vars <- unique(c(tabla_m1$Variable, tabla_m2$Variable, tabla_m3$Variable))

# Crear estructura base para la comparación
comparativa <- data.frame(Variable = todas_vars, stringsAsFactors = FALSE)

# Emparejar los OR de cada modelo a la variable correspondiente usando match()
comparativa$OR_M1 <- tabla_m1$OR[match(comparativa$Variable, tabla_m1$Variable)]
comparativa$OR_M2 <- tabla_m2$OR[match(comparativa$Variable, tabla_m2$Variable)]
comparativa$OR_M3 <- tabla_m3$OR[match(comparativa$Variable, tabla_m3$Variable)]

# Contar en cuántos modelos aparece cada variable para filtrar
comparativa$Conteo <- rowSums(!is.na(comparativa[, c("OR_M1", "OR_M2", "OR_M3")]))

# Filtrar solo las que aparecen en más de un modelo y excluir la columna auxiliar
comparativa_final <- comparativa[comparativa$Conteo > 1, c("Variable", "OR_M1", "OR_M2", "OR_M3")]

# Imprimir tabla comparativa
cat("\n--- RESUMEN COMPARATIVO: OR EN MÚLTIPLES MODELOS ---\n")
print(comparativa_final, row.names = FALSE)
