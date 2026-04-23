# 8. Inferencia del árbol de especies con ASTRAL

[← Volver al README](../README.md)

## 8.1 Preparar el conjunto de árboles génicos

```bash
Orchidaceae_allSub_63g_MLtree_20LBScoll.tre
```

## 8.2 Ejecutar ASTRAL

```bash
mkdir -p ASTRAL
java -jar astral.jar -i Orchidaceae_allSub_63g_MLtree_20LBScoll.tre -o ASTRAL/Orchidaceae_allSub_63g_ASTRAL.tre
```

## 8.3 Salida esperada

```bash
ASTRAL/Orchidaceae_allSub_63g_ASTRAL.tre
```
