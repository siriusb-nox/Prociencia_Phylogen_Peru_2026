# C. Alineamiento por locus y filtrado

[← Volver al README](../README.md)

El alineamiento múltiple de secuencias es el paso en el que se organizan secuencias homólogas de distintas especies para que sus posiciones sean comparables entre sí. La idea central es que cada columna del alineamiento represente una hipótesis de homología posicional, es decir, que los nucleótidos o aminoácidos en la misma columna descienden de una misma posición ancestral. Para lograr esto, los programas de alineamiento introducen gaps cuando ha habido inserciones o deleciones, de forma que las regiones comparables queden alineadas correctamente. El alineamiento es un paso crucial porque constituye la base de la inferencia filogenética: si el alineamiento es incorrecto, el árbol resultante también puede serlo. 

### Alineamiento de loci con `MAFFT`

Ejecuta los alineamientos por locus en paralelo:

```bash
ls *.FNA | parallel 'mafft --maxiterate 1000 --globalpair --adjustdirectionaccurately --thread 2 {} > ../alingSeqs_3xCov/{.}.fasta'
```

Este comando:

- recorre todos los loci recuperados con extensión `.FNA`
- alinea cada locus de forma independiente
- usa una estrategia de alineamiento global iterativa
- ajusta la orientación de las secuencias si es necesario

### 6.2 Filtrado de sitios mal alineados con `TAPER`

Crea primero la carpeta de salida:

```bash
mkdir -p Taper_filt
```

Después ejecuta:

```bash
for F in *.fasta; do julia correction_multi.jl -m N -a N -c 1 $F > Taper_filt/$F; done
```

Esto filtra sitios potencialmente mal alineados de cada alineamiento.

### 6.3 Eliminación de archivos vacíos

Tras el alineamiento o el filtrado, elimina archivos vacíos con:

```bash
find . -type f -empty
find . -type f -empty -delete
```

El primer comando muestra los archivos vacíos. El segundo los elimina.

### 6.4 Eliminación de loci con menos de 3 secuencias

Los loci con menos de tres secuencias no son adecuados para la inferencia de árboles génicos en este tutorial.

Este filtrado se ejecutó usando un script de shell ejecutable:

```bash
./filter_bySeqnum_v0.04102026.sh
```

Una versión simple de ese script es:

```bash
#!/usr/bin/env bash
set -euo pipefail

for f in *.fasta; do
    n=$(grep -c '^>' "$f")
    if [ "$n" -lt 3 ]; then
        rm "$f"
    fi
done
```

