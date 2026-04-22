"""
CONSTRUCCIÓN DE LA VARIABLE RESPUESTA Y

Este script toma la base inicial, y:

    (1) Clasifica cada fila como CUMPLIMIENTO / INCUMPLIMIENTO en DQO
        usando el corte de la escala de calidad de agua de CONAGUA

    (2) Construye la variable respuesta y (binaria: rojo = 1, no rojo = 0)
        que depende únicamente de DQO:
            - DQO incumple -> ROJO (y = 1)
            - DQO cumple   -> NO ROJO (y = 0)

SST y E-coli NO se usan para construir y; quedan disponibles como regresoras.
"""

import numpy as np
import pandas as pd

PATH_IN  = "/Users/danielarenee/Desktop/Water_Quality_Logit/base_filtrada.csv"
PATH_OUT = "base_con_y.csv"

# Umbral de cumplimiento
LIMITE_DQO = 40.0    # mg/L   (DQO ≤ 40  → cumple)

# Carga de la base
df = pd.read_csv(PATH_IN)
print(f"  Dimensiones: {df.shape}")

# Clasificación por parámetro
# Creamos una columna DQO_ESTADO con tres valores posibles:
# "CUMPLE" si valor ≤ límite, "INCUMPLE" si valor > límite, NaN si el valor original es faltante

def clasificar(serie, limite):
    resultado = pd.Series(np.nan, index=serie.index, dtype=object)
    resultado[serie <= limite] = "CUMPLE"
    resultado[serie >  limite] = "INCUMPLE"
    return resultado

df["DQO_ESTADO"] = clasificar(df["DQO_TOT"], LIMITE_DQO)

print(f"\n  DQO_ESTADO:")
for val, n in df["DQO_ESTADO"].value_counts(dropna=False).items():
    print(f"    {str(val):12s} {n:>5d}")

# Construcción de variable respuesta binaria
# y = 1  si DQO_ESTADO == "INCUMPLE"  (ROJO)
# y = 0  si DQO_ESTADO == "CUMPLE"    (NO ROJO)
# y = NaN si DQO_ESTADO es NaN        (no podemos decidir)

df["y"] = np.where(df["DQO_ESTADO"].isna(),
                   np.nan,
                   (df["DQO_ESTADO"] == "INCUMPLE").astype(float))

print("\n  Distribución de y:")
for val, n in df["y"].value_counts(dropna=False).items():
    pct = 100 * n / len(df)
    print(f"    {str(val):10s} {n:>5d}  ({pct:5.2f}%)")

n_validas = df["y"].notna().sum()
if n_validas > 0:
    prevalencia = df.loc[df["y"].notna(), "y"].mean()
    print(f"\n  Prevalencia de rojo (y = 1): {prevalencia:.2%}")
    print(f"  Filas con y bien definida:   {n_validas}/{len(df)}")

# Exportar
df.to_csv(PATH_OUT, index=False, encoding="utf-8-sig")
print(f"\n  ✔ Base exportada a: {PATH_OUT}")