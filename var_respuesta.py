"""
CONSTRUCCIÓN DE LA VARIABLE RESPUESTA Y

Este script toma la base inicial, y:

    (1) Clasifica cada fila como CUMPLIMIENTO / INCUMPLIMIENTO en los 3
        parámetros indicadores (DQO, SST y E-coli) usando los cortes de la
        escala de calidad de agua de CONAGUA

    (2) Construye la variable SEMAFORO (ROJO / AMARILLO / VERDE):
            - DQO incumple -> ROJO
            - DQO cumple, SST o EC incumplen -> AMARILLO
            - los tres cumplen -> VERDE

    (3) Construye la variable respuesta y (binaria: rojo = 1, no rojo = 0).
"""

import numpy as np
import pandas as pd

PATH_IN  = "/Users/danielarenee/Desktop/Water_Quality_Logit/base_filtrada.csv"
PATH_OUT = "base_con_y.csv"

# Umbrales de cumplimiento

LIMITE_DQO    =  40.0    # mg/L         (DQO ≤ 40  → cumple)
LIMITE_SST    = 150.0    # mg/L         (SST ≤ 150 → cumple)
LIMITE_E_COLI = 850.0    # NMP/100 mL   (EC  ≤ 850 → cumple)

# Carga de la base

df = pd.read_csv(PATH_IN)
print(f"  Dimensiones: {df.shape}")

# Clasificación por parámetro
# Para cada parámetro creamos una columna *_ESTADO con tres valores posibles:
# "CUMPLE" si valor ≤ límite, "INCUMPLE" si valor > límite, NaN si el valor original es faltante

def clasificar(serie, limite):
    resultado = pd.Series(np.nan, index=serie.index, dtype=object)
    resultado[serie <= limite] = "CUMPLE"
    resultado[serie >  limite] = "INCUMPLE"
    return resultado

df["DQO_ESTADO"] = clasificar(df["DQO_TOT"], LIMITE_DQO)
df["SST_ESTADO"] = clasificar(df["SST"],     LIMITE_SST)
df["EC_ESTADO"]  = clasificar(df["E_COLI"],  LIMITE_E_COLI)

for col in ["DQO_ESTADO", "SST_ESTADO", "EC_ESTADO"]:
    print(f"\n  {col}:")
    for val, n in df[col].value_counts(dropna=False).items():
        print(f"    {str(val):12s} {n:>5d}")

# Construcción del semáforo:

#     DQO incumple? ─── SÍ ─> ROJO
#          │ NO
#          v
#     DQO es NaN?  ─── SÍ ─> NaN  (no podemos decidir)
#          │ NO
#          v
#     SST o EC incumple? ─── SÍ ─> AMARILLO
#          │ NO
#          v
#     SST y EC cumplen? ─── SÍ ─> VERDE
#          │ NO (alguno es NaN, pero ninguno incumple)
#          v
#     NaN  (indeterminado entre amarillo y verde)

def semaforo(fila):
    dqo = fila["DQO_ESTADO"]
    sst = fila["SST_ESTADO"]
    ec  = fila["EC_ESTADO"]

    if dqo == "INCUMPLE": # Regla 1: DQO manda
        return "ROJO"
    if pd.isna(dqo): # Regla 2: si DQO es NaN, no podemos bajar a amarillo ni verde
        return np.nan
    # A partir de aquí, DQO cumple.
    if sst == "INCUMPLE" or ec == "INCUMPLE": # Regla 3: cualquier incumplimiento de SST o EC lo vuelve amarillo
        return "AMARILLO"
    if sst == "CUMPLE" and ec == "CUMPLE": # Regla 4: verde requiere que los tres cumplan (ninguno NaN)
        return "VERDE"
    # Caso residual: DQO cumple, ninguno incumple, pero SST o EC es NaN.
    # No se puede distinguir verde de amarillo → NaN.
    return np.nan

df["SEMAFORO"] = df.apply(semaforo, axis=1)

print("\n  Distribución del semáforo:")
for val, n in df["SEMAFORO"].value_counts(dropna=False).items():
    pct = 100 * n / len(df)
    print(f"    {str(val):10s} {n:>5d}  ({pct:5.2f}%)")

# Construcción de variable respuesta binaria
# y = 1  si SEMAFORO == "ROJO"
# y = 0  si SEMAFORO == "AMARILLO" o "VERDE"
# y = NaN si SEMAFORO es NaN

df["y"] = np.where(df["SEMAFORO"].isna(),
                   np.nan,
                   (df["SEMAFORO"] == "ROJO").astype(float))

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
