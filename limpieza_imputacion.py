"""
LIMPIEZA DE FILAS CON MUCHOS FALTANTES E IMPUTACIÓN MICE

Este script toma la base con binarias (salida de 03) y:

    (1) Elimina filas con y indefinida
    (2) Separa metadatos de variables regresoras candidatas
    (3) Excluye variables redundantes o colineales con y
    (4) Convierte TIPO CUERPO DE AGUA a dummy (LÓTICO = 1, LÉNTICO = 0).
    (5) Elimina filas con > UMBRAL_FILA de faltantes sobre las candidatas
    (6) Aplica MICE (Multiple Imputation by Chained Equations) para
        imputar los NaN restantes
    (7) Redondea las variables binarias imputadas a {0, 1}.

MICE funciona ajustando una regresión para cada variable con NaN usando
las demás como predictoras, e iterando hasta convergencia. Esto preserva
las relaciones entre las regresoras.
"""

import numpy as np
import pandas as pd

# IterativeImputer está aún en "experimental" en sklearn;
# hay que importarlo de forma especial.
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge

PATH_IN  = "/Users/danielarenee/Desktop/Water_Quality_Logit/base_con_binarias.csv"
PATH_OUT = "base_modelado.csv"

UMBRAL_FILA  = 0.50      # eliminar filas con > 50% de NaN en candidatas
MAX_ITER     = 10        # iteraciones de MICE
RANDOM_STATE = 42        # para reproducibilidad

# Metadatos: no entran al modelo, pero se conservan en la base de salida para trazabilidad
METADATOS = [
    "CLAVE SITIO", "CLAVE DE MONITOREO", "NOMBRE DEL SITIO",
    "FECHA REALIZACIÓN", "Año",
    "DQO_ESTADO", "SST_ESTADO", "EC_ESTADO", "SEMAFORO",
]

# Variables EXCLUIDAS del pool de candidatas a regresoras:
EXCLUIDAS = [
    "DQO_TOT", "DQO_SOL",              # construyen y
    "DBO_TOT", "DBO_SOL",              # casi colineales con DQO por construcción
    "COT", "COT_SOL",                  # ya recategorizadas como COT_NOM
    "TOX_V_15_UT", "TOX_V_5_UT",       # ya recategorizada como TOX_NOM
    "PB_TOT", "AS_TOT",                # casi sin variabilidad
    "PB_NOM", "AS_NOM",                # sus binarias también salen
]

# Variables binarias (tienen valores solo en {0, 1}) que hay que redondear
# después de imputar, porque MICE las trata como continuas
BINARIAS = ["COT_NOM", "TOX_NOM", "TIPO_LOTICO"]

# Carga y eliminación de filas sin y
df = pd.read_csv(PATH_IN)
print(f"  Base inicial: {df.shape}")

df = df[df["y"].notna()].reset_index(drop=True)
print(f"  Tras eliminar filas sin y: {df.shape}")

# Destokenizar COLOR_VER: "<5" → 2.5 porque se me fue jeje
df["COLOR_VER"] = (df["COLOR_VER"].astype(str).str.strip()
                   .replace("<5", "2.5")
                   .replace({"nan": None, "NaN": None, "": None}))
df["COLOR_VER"] = pd.to_numeric(df["COLOR_VER"], errors="coerce")

# dummy para cuerpo de agua
# LÓTICO = 1, LÉNTICO = 0
df["TIPO_LOTICO"] = (df["TIPO CUERPO DE AGUA"] == "LÓTICO").astype(int)
print(f"  TIPO_LOTICO = 1 (lóticos): {(df['TIPO_LOTICO']==1).sum()}")
print(f"  TIPO_LOTICO = 0 (lénticos): {(df['TIPO_LOTICO']==0).sum()}")

# separación de columnas
cols_excluir = set(METADATOS) | set(EXCLUIDAS) | {"y", "TIPO CUERPO DE AGUA"}
candidatas   = [c for c in df.columns if c not in cols_excluir]

print(f"  Candidatas a regresoras: {len(candidatas)}")
for c in candidatas:
    print(f"    - {c}")

# Eliminación de filas con umbral
# Calculamos el % de NaN por fila considerando solo las candidatas
pct_na_fila = df[candidatas].isna().mean(axis=1)
mask_conservar = pct_na_fila <= UMBRAL_FILA
df = df.loc[mask_conservar].reset_index(drop=True)

print(f"\n  Filas eliminadas: {(~mask_conservar).sum()}")
print(f"  Filas conservadas: {mask_conservar.sum()}")

# Reporte de NaNs
na_por_col = df[candidatas].isna().sum().sort_values(ascending=False)
total_na   = na_por_col.sum()
celdas     = len(df) * len(candidatas)
print(f"  Total de celdas a imputar: {total_na:,} "
      f"({100*total_na/celdas:.2f}% de la matriz)")
print(f"\n  NaN por columna (solo las que tienen):")
for c, n in na_por_col[na_por_col > 0].items():
    print(f"    {c:25s} {n:>5d}  ({100*n/len(df):5.2f}%)")


# Aplicación de MICE
# Usamos BayesianRidge como estimador base de cada regresión encadenada (opción por defecto de sklearn)
imputer = IterativeImputer(
    estimator     = BayesianRidge(),
    max_iter      = MAX_ITER,
    random_state  = RANDOM_STATE,
    verbose       = 0,
)

X_original = df[candidatas].copy()
X_imputado = imputer.fit_transform(X_original)
X_imputado = pd.DataFrame(X_imputado, columns=candidatas, index=df.index)

# Redondear binarias imputadas a {0, 1}
for b in BINARIAS:
    if b in X_imputado.columns:
        X_imputado[b] = (X_imputado[b] >= 0.5).astype(int)

# Reemplazar las candidatas imputadas en df
for c in candidatas:
    df[c] = X_imputado[c]

# Exportar
cols_finales = METADATOS + ["TIPO CUERPO DE AGUA"] + candidatas + ["y"]
cols_finales = [c for c in cols_finales if c in df.columns]
df_final = df[cols_finales].copy()

df_final.to_csv(PATH_OUT, index=False, encoding="utf-8-sig")
print("Base final...")
print(f"  Dimensiones: {df_final.shape}")
print(f"  Regresoras candidatas: {len(candidatas)}")
print(f"  Prevalencia de y = 1: {df_final['y'].mean():.2%}")
print(f"\n  ✔ Base exportada a: {PATH_OUT}")