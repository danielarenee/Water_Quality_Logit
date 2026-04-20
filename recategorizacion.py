"""
CATEGORIZACIÓN BINARIA DE PARÁMETROS SEGÚN LA NOM-001-SEMARNAT-2021

Este script crea versiones binarias de 4 parámetros numéricos,
según los Límites Permisibles (Promedio Diario) de la NOM-001-SEMARNAT-2021

    - Carbono Orgánico Total (COT)
    - Toxicidad Vibrio fischeri 15 min (TOX_V_15_UT)
    - Plomo total (PB_TOT)
    - Arsénico total (AS_TOT)

Los límites de COT y As dependen del tipo de cuerpo de agua:
    - LÓTICO  (ríos, arroyos, canales, drenes)
    - LÉNTICO (embalses, lagos y lagunas)

Una fila INCUMPLE la norma en un parámetro si el valor SUPERA el límite
correspondiente a su tipo. Para cada parámetro se crea una variable binaria:

    *_NOM : 0 = cumple, 1 = incumple, NaN si no clasificable

Estas variables binarias son las que se usarán como regresoras categóricas
en el modelo logit.
"""

import numpy as np
import pandas as pd

PATH_IN  = "/Users/danielarenee/Desktop/Water_Quality_Logit/base_con_y.csv"
PATH_OUT = "base_con_binarias.csv"

# Definición de las variables a categorizar.
#   clave: nombre de la columna binaria de salida
#   "origen" nombre de la columna numérica fuente en la base
#   "LÓTICO": límite de Promedio Diario para cuerpos lóticos
#   "LÉNTICO": límite de Promedio Diario para cuerpos lénticos
#
# Para agregar una variable nueva, añadir una entrada con el mismo formato

LIMITES_NOM = {
    "COT_NOM": {"origen": "COT",         "LÓTICO": 45.0, "LÉNTICO": 30.0 },
    "TOX_NOM": {"origen": "TOX_V_15_UT", "LÓTICO":  2.0, "LÉNTICO":  2.0 },
    "PB_NOM":  {"origen": "PB_TOT",      "LÓTICO":  0.3, "LÉNTICO":  0.3 },
    "AS_NOM":  {"origen": "AS_TOT",      "LÓTICO":  0.3, "LÉNTICO":  0.15},
}

# Carga
df = pd.read_csv(PATH_IN)
print(f"  Dimensiones: {df.shape}")

# Categorizar variable por variable
for nombre_salida, config in LIMITES_NOM.items():
    col_origen  = config["origen"]
    lim_lotico  = config["LÓTICO"]
    lim_lentico = config["LÉNTICO"]

    # serie con el límite aplicable a cada fila según su tipo
    limite_fila = df["TIPO CUERPO DE AGUA"].map(
        {"LÓTICO": lim_lotico, "LÉNTICO": lim_lentico}
    )

    # máscara de filas clasificables
    serie  = df[col_origen]
    valido = serie.notna() & limite_fila.notna()

    # Versión binaria (0 = cumple, 1 = incumple, NaN si no clasificable)
    binaria = pd.Series(np.nan, index=df.index, dtype=float)
    binaria[valido & (serie <= limite_fila)] = 0.0
    binaria[valido & (serie >  limite_fila)] = 1.0
    df[nombre_salida] = binaria

    # Reporte
    n_cumple   = int((binaria == 0).sum())
    n_incumple = int((binaria == 1).sum())
    n_na       = int(binaria.isna().sum())
    n_total    = len(df)
    print(f"\n  {col_origen}  →  {nombre_salida}  "
          f"(límites: lótico={lim_lotico}, léntico={lim_lentico})")
    print(f"    CUMPLE   (0) {n_cumple:>5d}  ({100*n_cumple/n_total:5.2f}%)")
    print(f"    INCUMPLE (1) {n_incumple:>5d}  ({100*n_incumple/n_total:5.2f}%)")
    print(f"    NaN          {n_na:>5d}  ({100*n_na/n_total:5.2f}%)")

# Resumen
nombres_binarias = list(LIMITES_NOM.keys())
resumen = pd.DataFrame({
    "n_cumple":    [(df[b] == 0).sum()                 for b in nombres_binarias],
    "n_incumple":  [(df[b] == 1).sum()                 for b in nombres_binarias],
    "n_faltante":  [df[b].isna().sum()                 for b in nombres_binarias],
    "pct_incumple":[100 * (df[b] == 1).mean().round(4) for b in nombres_binarias],
}, index=nombres_binarias)
print(resumen.to_string())

# Exportar
df.to_csv(PATH_OUT, index=False, encoding="utf-8-sig")
print(f"\n  ✔ Base exportada a: {PATH_OUT}")
print(f"  Dimensiones finales: {df.shape}")