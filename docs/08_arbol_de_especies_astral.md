# Inferencia del árbol de especies con ASTRAL

[← Volver al README](../README.md)

ASTRAL es un método para inferir un árbol de especies a partir de un conjunto de árboles génicos. En lugar de concatenar todas las secuencias en una sola supermatriz, ASTRAL resume la información contenida en muchos gene trees y busca el árbol de especies que sea más congruente con ellos. Lo hace comparando cuartetos —subárboles de cuatro taxones— y seleccionando la topología que maximiza el acuerdo entre el árbol de especies y los árboles génicos de entrada. ASTRAL fue diseñado dentro del marco del multispecies coalescent, por lo que es especialmente útil cuando existe discordancia entre genes debida a procesos como incomplete lineage sorting (ILS).  

En términos simples, ASTRAL toma muchos árboles de genes, los resume y produce una hipótesis de árbol de especies que refleja la señal compartida entre loci, incluso cuando los árboles génicos no son idénticos entre sí. También puede calcular soporte local de ramas a partir de la información de cuartetos, lo que ayuda a interpretar qué relaciones están mejor respaldadas por los datos.


Crea primero una carpeta para guardar la salida:

```bash

mkdir -p ASTRAL

```

Después ejecuta ASTRAL:

```bash

java -jar astral.jar -i Orchidaceae_allSub_63g_MLtree_20LBScoll.tre -o ASTRAL/Orchidaceae_allSub_63g_ASTRAL.tre

```

Este comando:

- usa como entrada el conjunto de árboles génicos filtrados

- infiere un **árbol de especies** bajo el marco del **multispecies coalescent**

- escribe el árbol final en formato Newick

### Interpretación general del resultado

El árbol resultante:

- representa una **hipótesis de relaciones entre especies**

- resume la señal compartida entre múltiples árboles génicos

- puede ayudar a identificar conflicto genealógico y discordancia entre loci

### Salida esperada

El archivo principal será:

```bash

ASTRAL/Orchidaceae_allSub_63g_ASTRAL.tre

```

Este árbol puede visualizarse posteriormente en programas como **FigTree**, **iTOL** o **Dendroscope**.

### Versión opcional con archivo de log

Si quieres guardar también el log de ejecución, puedes correr:

```bash

java -jar astral.jar \\

  -i Orchidaceae_allSub_63g_MLtree_20LBScoll.tre \\

  -o ASTRAL/Orchidaceae_allSub_63g_ASTRAL.tre \\

  > ASTRAL/astral.log

```

### ¿Por qué usar ASTRAL?

ASTRAL es especialmente útil cuando:

- distintos genes apoyan historias evolutivas diferentes

- existe discordancia entre árboles génicos

- se desea inferir un árbol de especies a partir de múltiples loci independientes

En este contexto, ASTRAL complementa la inferencia de árboles génicos por **RAxML-NG** y permite pasar de una colección de gene trees a una hipótesis de árbol de especies.
