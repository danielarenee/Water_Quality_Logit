"""
SPLIT TRAIN/TEST, INSPECCIÓN DE DISTRIBUCIONES Y MATRIZ DE CORRELACIONES

Este script toma la base modelada (salida de 04) y:

    (1) Separa la base en 80% entrenamiento y 20% prueba, de forma
        estratificada por y (aka. mantiene la misma proporción ROJO/NO
        ROJO en ambos subconjuntos)
    (2) Sobre el 80% de entrenamiento:
        (a) calcula estadísticos descriptivos y sesgo de cada continua
        (b) calcula la matriz de correlaciones de Spearman
        (c) Identifica pares con |correlación| > 0.7 (candidatas a ser
            familias redundantes)
        (d) Grafica la matriz de correlaciones como heatmap

El 20% de prueba se guarda intacto y no se usa para ninguna selección.
Solo se usará al final para la matriz de confusión. Se hace desde ahorita
para evitar data leaking

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

PATH_IN         = "base_modelado.csv"
PATH_TRAIN      = "base_train.csv"
PATH_TEST       = "base_test.csv"
PATH_CORR       = "matriz_correlaciones_train.csv"
PATH_HEATMAP    = "heatmap_correlaciones.png"

TEST_SIZE       = 0.20
RANDOM_STATE    = 42 # para reproducibilidad

# umbrales para diagnóstico
UMBRAL_SESGO    = 1.0
UMBRAL_CORR     = 0.70

# metadatos y binarias
METADATOS = [
    "CLAVE SITIO", "CLAVE DE MONITOREO", "NOMBRE DEL SITIO",
    "FECHA REALIZACIÓN", "Año",
    "DQO_ESTADO", "SST_ESTADO", "EC_ESTADO", "SEMAFORO",
    "TIPO CUERPO DE AGUA",
]
BINARIAS = ["COT_NOM", "TOX_NOM", "TIPO_LOTICO"]

# Carga y split estratificado

df = pd.read_csv(PATH_IN)
print(f"  Base completa: {df.shape}")
print(f"  Prevalencia de y=1 (total): {df['y'].mean():.3f}")

train, test = train_test_split(
    df,
    test_size    = TEST_SIZE,
    stratify     = df["y"],
    random_state = RANDOM_STATE,
)
train = train.reset_index(drop=True)
test  = test.reset_index(drop=True)

print(f"\n  Train (80%): {train.shape}   prevalencia y=1: {train['y'].mean():.3f}")
print(f"  Test  (20%): {test.shape}    prevalencia y=1: {test['y'].mean():.3f}")

# Exportar test
train.to_csv(PATH_TRAIN, index=False, encoding="utf-8-sig")
test.to_csv(PATH_TEST,   index=False, encoding="utf-8-sig")
print(f"\n  ✔ Train exportado a: {PATH_TRAIN}")
print(f"  ✔ Test  exportado a: {PATH_TEST}  (CONGELADO, no tocar)")


# Identificar continuas y binarias
continuas = [c for c in train.columns
             if c not in METADATOS + BINARIAS + ["y"]
             and pd.api.types.is_numeric_dtype(train[c])]

print("Clasificación de columnas... ")
print(f"  Variables continuas:  {len(continuas)}")
print(f"  Variables binarias:   {len(BINARIAS)}")
print(f"  Metadatos:            {len(METADATOS)}")

# Ahora, en train... EDA
# el sesgo (skewness) mide la asimetría de la distribución. tomaremos >1 como sesgo fuerte

print("Estadísticos y sesgos...")
desc = pd.DataFrame({
    "min":       train[continuas].min().round(3),
    "mediana":   train[continuas].median().round(3),
    "media":     train[continuas].mean().round(3),
    "max":       train[continuas].max().round(3),
    "skewness":  train[continuas].skew().round(3),
})
desc["sesgo_fuerte"] = desc["skewness"].abs() > UMBRAL_SESGO
desc = desc.sort_values("skewness", ascending=False)

print(desc.to_string())

# En train... matriz de correlaciones de Spearman
# Spearman porque es robusta a outliers, captura relaciones no lineales,
# es invariante a transformaciones monótonas

print("Matriz de correlaciones...")
variables_corr = continuas + BINARIAS
corr = train[variables_corr].corr(method="spearman")

# buscamos los pares altamente correlacionados
print(f"\n  Pares con |correlación Spearman| > {UMBRAL_CORR}:")
pares_altos = []
for i in range(len(corr.columns)):
    for j in range(i + 1, len(corr.columns)):
        r = corr.iloc[i, j]
        if abs(r) > UMBRAL_CORR:
            pares_altos.append((corr.columns[i], corr.columns[j], r))

pares_altos.sort(key=lambda x: -abs(x[2]))
for v1, v2, r in pares_altos:
    print(f"    {v1:20s} <-> {v2:20s}  r = {r:+.3f}")

# exportar matriz
corr.round(3).to_csv(PATH_CORR, encoding="utf-8-sig")
print(f"\n  ✔ Matriz completa exportada a: {PATH_CORR}")

# graficar heatmap de la matriz de correlaciones
# rojo = correlación positiva fuerte, azul = negativa fuerte, blanco = cero

fig, ax = plt.subplots(figsize=(14, 12))
im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")

# etiquetas en ambos ejes
ax.set_xticks(np.arange(len(corr.columns)))
ax.set_yticks(np.arange(len(corr.columns)))
ax.set_xticklabels(corr.columns, rotation=90, fontsize=8)
ax.set_yticklabels(corr.columns, fontsize=8)

# escribir el valor numérico dentro de cada celda (solo si |r| > 0.5 para no saturar)
for i in range(len(corr.columns)):
    for j in range(len(corr.columns)):
        r = corr.values[i, j]
        if abs(r) > 0.5 and i != j:
            color = "white" if abs(r) > 0.75 else "black"
            ax.text(j, i, f"{r:.2f}", ha="center", va="center",
                    color=color, fontsize=6)

ax.set_title("Matriz de correlaciones de Spearman (train)", fontsize=13, pad=15)
fig.colorbar(im, ax=ax, shrink=0.8, label="Correlación")
fig.tight_layout()
fig.savefig(PATH_HEATMAP, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  ✔ Heatmap guardado en: {PATH_HEATMAP}")

# hacemos un screening preliminar viendo la correlación de cada regresora candidata con y

print("Correlaciones con y...")

corr_y = train[variables_corr + ["y"]].corr(method="spearman")["y"].drop("y")
corr_y = corr_y.abs().sort_values(ascending=False)

print(f"\n  Ordenadas por |correlación con y|:")
for v, r in corr_y.items():
    signo = "+" if train[[v, "y"]].corr(method="spearman").iloc[0, 1] >= 0 else "-"
    print(f"    {v:20s}  |r| = {r:.3f}  (signo: {signo})")