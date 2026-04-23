# ==============================================================================
# 01_seleccion_variables.R
# ==============================================================================
# SELECCIÓN DE VARIABLES VÍA VIF (EN R)
#
# Continuación del pipeline que veníamos trabajando en Python. Este script:
#
#   (1) Carga base_train.csv y base_test.csv (ya creados por Python).
#   (2) Verifica la matriz de correlaciones de Spearman (debe coincidir
#       con la que calculó Python).
#   (3) Grafica el heatmap de correlaciones.
#   (4) Elimina las 2 variables "gemelas perfectas" que decidimos sacar
#       a mano: SDT (r=1.00 con CONDUC_CAMPO) y OD_% (r=0.95 con OD_mg/L).
#   (5) Aplica VIF iterativo INFORMADO para eliminar multicolinealidad:
#       en cada iteración identifica todas las variables con VIF > umbral,
#       y de entre ellas elimina la que MENOS correlaciona con y (Spearman).
#       Así no eliminamos una variable solo por su redundancia si resulta
#       que es la que mejor predice y dentro de su familia.
#   (6) Exporta la base con el pool final de regresoras listas para
#       proponer los 3 modelos.
#
# Requisitos de paquetes:
#   install.packages(c("car", "corrplot", "dplyr"))
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. CONFIGURACIÓN
# ------------------------------------------------------------------------------
# Ajusta esta ruta a donde tengas los CSV
setwd("/Users/danielarenee/Desktop/Water_Quality_Logit")

library(car)       # para vif()
library(corrplot)  # para graficar la matriz de correlaciones
library(dplyr)

PATH_TRAIN   <- "base_train.csv"
PATH_TEST    <- "base_test.csv"
PATH_OUT     <- "base_train_final.csv"
PATH_OUT_TS  <- "base_test_final.csv"
PATH_HEATMAP <- "heatmap_correlaciones_R.png"

UMBRAL_VIF   <- 10   # VIF > 10 se considera multicolinealidad crítica

# Variables que NO entran al modelo:
#   - metadatos (identificadores, fechas, estados intermedios)
#   - SDT  (r = 1.00 con CONDUC_CAMPO; decisión manual)
#   - OD_% (r = 0.95 con OD_mg/L;      decisión manual)
METADATOS <- c("CLAVE SITIO", "CLAVE DE MONITOREO", "NOMBRE DEL SITIO",
               "FECHA REALIZACIÓN", "Año",
               "DQO_ESTADO", "SST_ESTADO", "EC_ESTADO", "SEMAFORO",
               "TIPO CUERPO DE AGUA")

ELIMINAR_MANUAL <- c("SDT", "OD_.")   # OD_% llega a R como OD_. por el %


# ------------------------------------------------------------------------------
# 1. CARGA DE DATOS
# ------------------------------------------------------------------------------
cat("=====================================================\n")
cat("PASO 1: CARGA DE DATOS\n")
cat("=====================================================\n")

train <- read.csv(PATH_TRAIN, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")
test  <- read.csv(PATH_TEST,  check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")

# R interpreta "%" en nombres de columna como carácter raro, renombramos
names(train)[names(train) == "OD_%"] <- "OD_pct"
names(test)[names(test)  == "OD_%"]  <- "OD_pct"
ELIMINAR_MANUAL <- c("SDT", "OD_pct")

cat(sprintf("  Train: %d filas x %d columnas\n", nrow(train), ncol(train)))
cat(sprintf("  Test:  %d filas x %d columnas\n", nrow(test),  ncol(test)))
cat(sprintf("  Prevalencia y=1 en train: %.3f\n", mean(train$y)))
cat(sprintf("  Prevalencia y=1 en test:  %.3f\n", mean(test$y)))


# ------------------------------------------------------------------------------
# 2. IDENTIFICACIÓN DEL POOL DE CANDIDATAS
# ------------------------------------------------------------------------------
cat("\n=====================================================\n")
cat("PASO 2: POOL INICIAL DE CANDIDATAS\n")
cat("=====================================================\n")

# candidatas = todas las numéricas menos metadatos, y, y las que eliminamos a mano
cols_excluir <- c(METADATOS, "y", ELIMINAR_MANUAL)
candidatas   <- setdiff(names(train), cols_excluir)
# solo numéricas
candidatas   <- candidatas[sapply(train[candidatas], is.numeric)]

cat(sprintf("  Candidatas iniciales: %d\n", length(candidatas)))
cat(sprintf("  Eliminadas manualmente: %s\n",
            paste(ELIMINAR_MANUAL, collapse = ", ")))


# ------------------------------------------------------------------------------
# 3. MATRIZ DE CORRELACIONES DE SPEARMAN (VERIFICACIÓN CON PYTHON)
# ------------------------------------------------------------------------------
# Debe coincidir en magnitudes con la que calculó Python. Pequeñas diferencias
# en el 3er decimal son normales (diferentes implementaciones de empates).
cat("\n=====================================================\n")
cat("PASO 3: MATRIZ DE CORRELACIONES\n")
cat("=====================================================\n")

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


# ------------------------------------------------------------------------------
# 4. HEATMAP
# ------------------------------------------------------------------------------
cat("\n=====================================================\n")
cat("PASO 4: HEATMAP DE CORRELACIONES\n")
cat("=====================================================\n")

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
         col       = colorRampPalette(c("#08306B", "white", "#67000D"))(200),
         title     = "Matriz de correlaciones de Spearman (train, R)",
         mar       = c(0, 0, 2, 0))
dev.off()
cat(sprintf("  ✔ Heatmap guardado en: %s\n", PATH_HEATMAP))


# ------------------------------------------------------------------------------
# 5. VIF ITERATIVO INFORMADO
# ------------------------------------------------------------------------------
# El VIF mide qué tanto se infla la varianza del coeficiente de una variable
# por su correlación con las demás. Se define como:
#       VIF_i = 1 / (1 - R²_i)
# donde R²_i es el R² de regresar x_i contra el resto de las regresoras.
#
# Importante: VIF NO considera a y. Mide redundancia entre X, no utilidad
# predictiva. Por eso, si solo elimináramos "la variable de mayor VIF",
# podríamos estar tirando justo la más útil de una familia redundante.
#
# Procedimiento INFORMADO (el que usamos):
#   - calcular VIF de todas las variables
#   - identificar todas las que tengan VIF > umbral
#   - DE ENTRE ESAS, eliminar la que menos correlacione con y (Spearman)
#   - recalcular, repetir hasta que todas tengan VIF <= umbral
#
# Así combinamos el requisito de bondad de ajuste (VIF <= 10, pedido por
# la rúbrica) con la consideración de que la variable eliminada debe ser
# la menos útil, no la más redundante.
#
# Nota técnica: vif() de car necesita un modelo (lm o glm). Como es un
# diagnóstico entre regresoras, usamos un lm auxiliar con y como
# respuesta. La selección NO depende del modelo, solo de las X entre sí.

cat("\n=====================================================\n")
cat("PASO 5: VIF ITERATIVO INFORMADO (umbral = 10)\n")
cat("=====================================================\n")

# funciones auxiliares
calcular_vif <- function(vars, data) {
  formula <- as.formula(paste("y ~", paste(paste0("`", vars, "`"),
                                           collapse = " + ")))
  modelo  <- lm(formula, data = data)
  vif(modelo)
}

# pre-cómputo: correlación de Spearman de cada candidata con y
# (no cambia en el proceso: y no se toca)
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
  
  # entre las culpables, elegir la de MENOR correlación con y
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

cat("\n-----------------------------------------------------\n")
cat("VIF FINAL DE LAS VARIABLES SUPERVIVIENTES\n")
cat("-----------------------------------------------------\n")
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


# ------------------------------------------------------------------------------
# 6. EXPORTACIÓN DEL POOL FINAL
# ------------------------------------------------------------------------------
cat("\n=====================================================\n")
cat("PASO 6: EXPORTACIÓN DEL POOL FINAL\n")
cat("=====================================================\n")

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

cat("\n=====================================================\n")
cat("FIN. Próximo paso: proponer los 3 modelos candidatos.\n")
cat("=====================================================\n")