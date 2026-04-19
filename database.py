"""
CONSTRUCCIÓN DE LA BASE DE DATOS

Este script realiza lo siguiente:
    (1) Carga las 3 pestañas del archivo .xlsb
    (2) Limpieza de strings en columnas categóricas (quita espacios invisibles)
    (3) Filtrado de filas: solo cuerpos de agua LÓTICOS y LÉNTICOS
        (descarta humedales, subterráneos y costeros)
    (4) Conversión de columnas tipo texto a numéricas, manejando valores
        censurados del tipo "<L" (debajo del límite de detección) como L/2
    (5) Eliminación de columnas con <70% de cobertura de datos
"""

import numpy as np
import pandas as pd

PATH_XLSB          = "/Users/danielarenee/Desktop/Water_Quality_Logit/datos_logit.xlsb"
PATH_OUT           = "base_filtrada.csv"
UMBRAL_COBERTURA   = 0.70 # 70% mínimo de valores no nulos

# Columnas que NO deben convertirse a numéricas
# (identificadores, fechas y campos textuales)
COLS_NO_NUMERICAS = {
    "CLAVE SITIO", "CLAVE DE MONITOREO", "NOMBRE DEL SITIO",
    "TIPO CUERPO DE AGUA", "FECHA REALIZACIÓN", "Año",
    "MAT_FLOTANTE", "COLOR_VER",
}

# 1. Carga de las pestañas

sitios = pd.read_excel(PATH_XLSB, sheet_name="Sitos",      engine="pyxlsb")
datos  = pd.read_excel(PATH_XLSB, sheet_name="Resultados", engine="pyxlsb")
simbo  = pd.read_excel(PATH_XLSB, sheet_name="Simbología", engine="pyxlsb")

print(f"  Sitos:      {sitios.shape}")
print(f"  Resultados: {datos.shape}")
print(f"  Simbología: {simbo.shape}")


# 2. Limpieza de strings en columnas categoricas

datos["TIPO CUERPO DE AGUA"]      = datos["TIPO CUERPO DE AGUA"].astype(str).str.strip()
sitios["TIPO DE CUERPO DE AGUA"]  = sitios["TIPO DE CUERPO DE AGUA"].astype(str).str.strip()
sitios["SUBTIPO CUERPO AGUA"]     = sitios["SUBTIPO CUERPO AGUA"].astype(str).str.strip()

print("  Distribución en Resultados ANTES del filtrado:")
for tipo, n in datos["TIPO CUERPO DE AGUA"].value_counts().items():
    print(f"    {tipo:25s} {n:>5d}")

# 3. Filtrado de lótico y léntico

df = datos[datos["TIPO CUERPO DE AGUA"].isin(["LÓTICO", "LÉNTICO"])].copy()
df = df.reset_index(drop=True)

print(f"  Observaciones: {len(datos)} → {len(df)}  "
      f"(descartadas: {len(datos) - len(df)})")
for tipo, n in df["TIPO CUERPO DE AGUA"].value_counts().items():
    print(f"    {tipo:25s} {n:>5d}")

# 4. Destokenización y corrección de NaNs manuales

def destokenizar(x):
    """Convierte valores tipo texto (p. ej. '<0.015') a flotante.
    - '<L' → L/2   (censurado a la izquierda: por debajo de LOD)
    - '>L' → L     (saturado: por encima del rango)
    - 'ND', 'NA', 'AUSENTE', ''  → NaN
    - Strings no parseables → NaN
    """
    if pd.isna(x):
        return np.nan
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.upper() in ("ND", "NA", "AUSENTE"):
        return np.nan
    try:
        if s.startswith("<"):
            return float(s[1:]) / 2.0
        if s.startswith(">"):
            return float(s[1:])
        return float(s)
    except ValueError:
        return np.nan

cols_a_convertir = [c for c in df.columns
                    if c not in COLS_NO_NUMERICAS and df[c].dtype == object]
print(f"  Columnas convertidas a numéricas: {len(cols_a_convertir)}")

for col in cols_a_convertir:
    df[col] = df[col].apply(destokenizar)


# 5. Filtro de columnas por cobertura (>=70%)
# Esto para que los metodos de imputación que eligamos para los NaNs sean fiables

cobertura = df.notna().mean()
cols_a_conservar = cobertura[cobertura >= UMBRAL_COBERTURA].index.tolist()
cols_a_eliminar  = cobertura[cobertura <  UMBRAL_COBERTURA].index.tolist()

print(f"  Columnas originales: {df.shape[1]}")
print(f"  Columnas conservadas: {len(cols_a_conservar)}")
print(f"  Columnas eliminadas:  {len(cols_a_eliminar)}")

df_final = df[cols_a_conservar].copy()


# Reporte y exportación inicial de la base
print("REPORTE FINAL")
print(f"  Dimensiones de la base filtrada: {df_final.shape}")
print(f"\n  Cobertura de cada columna conservada (ordenada):")
cob_final = (df_final.notna().mean().mul(100).round(2)
             .sort_values(ascending=False))
for c, v in cob_final.items():
    print(f"    {c:25s} {v:6.2f}%")

df_final.to_csv(PATH_OUT, index=False, encoding="utf-8-sig")
print(f"\n  ✔ Base exportada a: {PATH_OUT}")
