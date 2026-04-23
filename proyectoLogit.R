# ==============================================================================
# PROYECTO DE REGRESIÓN LOGÍSTICA
# DANIELA RENÉE Y MARIO SÁNCHEZ
# ECONOMETRÍA 1
# ==============================================================================

setwd("/Users/danielarenee/Desktop/Water_Quality_Logit")

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

# aqui agregar el base 