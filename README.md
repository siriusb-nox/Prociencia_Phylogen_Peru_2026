# Taller de "Principios en Filogenómica" - 15-19 Abril (Madre de Dios) y 22-24 Abril (UNTRM) 
## Instructors: [Oscar Alejandro Perez-Escobar](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.70484) [(Royal Botanic Gardens, Kew)](https://scholar.google.com/citations?user=tSzyp6QAAAAJ&hl=es) & [Alexander Damian](https://cameron.botany.wisc.edu/staff/damian-l-alexander/) [(University of Wisconsin-Madison)](https://scholar.google.co.uk/citations?user=RwjJKcwAAAAJ&hl=th)
##### Este tayer esta apoyado por [Prociencia](https://prociencia.gob.pe/) & la [Swiss Orchid Foundation](https://novagenetix.ch/swiss-orchid-foundation/)
## Introduccion
Este repositorio contiene un tutorial introductorio para una práctica de filogenómica ejecutada en una máquina virtual de Ubuntu sobre portátiles con Windows. El tutorial cubre un flujo de trabajo mínimo para descargar lecturas Illumina, limpiar los datos, ensamblar loci in-silico, alinear genes, filtrar loci y sitios problemáticos, inferir árboles génicos de máxima verosimilitud y arboles de especies usando un metodo de "Multispecies Coalescent" (MSC).

En este taller, emplearemos una pequeña muestra de datos de la filogenia de orquideas de Perez-Escobar et al. (2024, [_New. Phyt._](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.19580)), que representa las cinco subfamilias mas un grupo externo (_Dioscorea_). El objetivo del taller es lograr ensamblar una filogenia a partir de sets de datos illumina los cuales primero deben ser obtenidos de la base de datos SRA del [NCBI](https://www.ncbi.nlm.nih.gov/) asi como tambien aprender a diseñar proyectos de investigacion en filogenomica. 

<p align="center">
  <img src="https://github.com/siriusb-nox/Prociencia_Phylogen_Peru_2026/blob/main/IMG/Orchid_ToL3.0_v5_OP_distproof_LR_96dpi_MQ_online_EN.jpg" alt="Orchid Tree of Life v2.0"/>
</p>


---

## 1. Entorno computacional del taller

La práctica fue diseñada para ejecutarse en una máquina virtual de Ubuntu en VirtualBox. Alternativamente, los ejercicios del taller se podran implementar en un ambiente Linux, asumiendo que todas las dependencias han sido instaladas.

### 1.1 Información de la máquina virtual

- **Sistema operativo invitado:** Ubuntu
- **Software de virtualización:** VirtualBox
- **Usuario:** `vboxuser`

> **La contraseña de cada sesion Ubuntu sera entregada el dia de la practica**.

### 1.2 Programas requeridos
Para ejecutar este taller, se requieren los siguientes programas:

- [VirtualBox](https://www.virtualbox.org/wiki/Downloads)
- Ubuntu (máquina virtual)
- Bash
- [SRA Toolkit: `prefetch` y `fasterq-dump`](https://github.com/ncbi/sra-tools)
- [fastp](https://github.com/opengene/fastp)
- [conda](https://docs.conda.io/)
- [HybPiper](https://github.com/mossmatters/HybPiper)
- [BWA](https://github.com/lh3/BWA)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- [Julia](https://julialang.org/)
- [TAPER](https://github.com/chaoszhang/TAPER)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- [Newick Utilities](https://github.com/tjunier/newick_utils)
- utilidades estándar de UNIX: `grep`, `find`, `rename`, `cat`, `mv`, `ls`, `gzip`

---

## 2. Introducción muy breve a Bash

Esta sección está adaptada del tutorial de Bash usado en el repositorio del taller de ONT y simplificada para esta práctica de filogenómica.

### 2.1 ¿Qué es Bash?

Bash es la shell usada en muchos sistemas operativos de tipo UNIX, como Ubuntu y Linux. En bioinformática, la shell es especialmente útil porque muchas herramientas están diseñadas para ejecutarse desde la terminal y porque Bash facilita el trabajo con archivos de texto, secuencias y directorios.

Un comando típico tiene esta estructura:

```bash
programa opcion1 opcion2 entrada > salida
```

Por ejemplo:

```bash
ls
ls -lrt
```

El primer comando lista archivos. El segundo los muestra con más detalle y los ordena por fecha.

### 2.2 Comandos básicos

#### `pwd`

Muestra el directorio actual:

```bash
pwd
```

#### `ls`

Lista los archivos del directorio actual:

```bash
ls
ls -l
```

#### `cd`

Cambia de directorio:

```bash
cd
cd ..
cd /home/vboxuser/Documents
```

#### `mkdir`

Crea directorios:

```bash
mkdir nueva_carpeta
mkdir carpeta1 carpeta2 carpeta3
mkdir -p ~/Documents/Phylo_work_2026
```

#### `cat`

Muestra o concatena archivos de texto:

```bash
cat archivo.txt
cat archivo1.txt archivo2.txt > combinado.txt
```

#### `less`

Permite inspeccionar archivos de texto página por página:

```bash
less archivo.txt
```

#### `grep`

Busca patrones de texto:

```bash
grep '^>' archivo.fasta
grep -c '^>' archivo.fasta
```

#### `man`

Muestra la página de ayuda de un comando:

```bash
man ls
man grep
```

### 2.3 Ideas útiles de Bash para este taller

#### Bucles

Un bucle repite un comando para muchos archivos o muestras:

```bash
for F in *.fasta; do echo "$F"; done
```

#### Redirección de salida

El símbolo `>` envía la salida a un archivo:

```bash
cat entrada.txt > salida.txt
```

#### Tuberías

Las tuberías envían la salida de un programa a otro:

```bash
ls *.FNA | wc -l
```

#### Ejecución en paralelo

GNU Parallel permite correr muchos trabajos al mismo tiempo:

```bash
ls *.FNA | parallel 'echo {}'
```

---

## 3. Estructura sugerida de directorios de trabajo

Una estructura ordenada de carpetas hace que la práctica sea mucho más fácil de seguir. 

```bash
mkdir -p /media/foldername/Phylo_work_2026_Peru/{raw_reads,trimmed,HybPiper,alignSeqs_3xCov,Taper_filt,raxml_gene_trees}
```

También puedes crear carpetas separadas para scripts, logs y archivos de referencia.

---

## 4. Descarga y limpieza de datos Illumina

### 4.1 Descarga de lecturas crudas desde NCBI SRA

Prepara un archivo de texto plano con un acceso por línea, por ejemplo `accessions.txt`:

```text
ERR5006156
ERR7622624
ERR7622662
SRR26876932
SRR26876810
SRR26876866
SRR6127594
SRR6127593
ERR7599286
ERR4180009
```

Guarda el siguiente script como `download_sra_fastq.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

ACC_FILE="${1:-accessions.txt}"
OUTDIR="${2:-$PWD/sra_fastq}"
THREADS="${3:-4}"

TMPDIR="${OUTDIR}/tmp"
SRA_DIR="${OUTDIR}/sra_cache"
FASTQ_DIR="${OUTDIR}/fastq"

mkdir -p "${TMPDIR}" "${SRA_DIR}" "${FASTQ_DIR}"

while IFS= read -r ACC || [[ -n "$ACC" ]]; do
    [[ -z "$ACC" ]] && continue

    echo "Processing ${ACC} ..."

    prefetch "${ACC}" --output-directory "${SRA_DIR}"

    fasterq-dump "${SRA_DIR}/${ACC}" \
        --split-files \
        --threads "${THREADS}" \
        --temp "${TMPDIR}" \
        --outdir "${FASTQ_DIR}"

    gzip -f "${FASTQ_DIR}/${ACC}"*.fastq

    echo "Finished ${ACC}"
done < "${ACC_FILE}"

echo "All downloads completed."
echo "Compressed FASTQ files are in: ${FASTQ_DIR}"
```

Ejecútalo así:

```bash
chmod +x download_sra_fastq.sh
./download_sra_fastq.sh accessions.txt ~/Documents/Phylo_work_2026/raw_reads 4
```

### 4.2 Renombrar archivos FASTQ con nombres de especies y códigos de acceso

Supón que tienes un archivo tabulado de dos columnas llamado `accession_to_species.tsv`:

```text
ERR5006156	Apostasia_nuda
ERR7622624	Aa_maderoi
ERR7622662	Anacamptis_pyramidalis
SRR26876932	Cycnoches_amparoanum
SRR26876810	Coryanthes_hunteriana
SRR26876866	Vanilla_hartii
SRR6127594	Phragmipedium_lindleyanum
SRR6127593	Paphiopedilum_callosum
ERR7599286	Erythrorchis_cassythoides
ERR4180009	Dioscorea_caucasica
```

Puedes renombrar los FASTQ descargados con un script de Bash como este:

```bash
#!/usr/bin/env bash
set -euo pipefail

MAPFILE="accession_to_species.tsv"
declare -A species_map

while IFS=$'\t' read -r acc species; do
    [[ -z "${acc}" ]] && continue
    species_map["$acc"]="$species"
done < "$MAPFILE"

for f in *_*.fastq.gz; do
    [[ -e "$f" ]] || continue

    acc="${f%%_*}"
    rest="${f#*_}"
    pair="${rest%%.*}"

    if [[ -n "${species_map[$acc]:-}" ]]; then
        new="${species_map[$acc]}_${acc}_${pair}.fastq.gz"
        mv -- "$f" "$new"
    fi
done
```

### 4.3 Recorte de lecturas con `fastp`

El conjunto final de muestras usado en la práctica tras el recorte fue:

```text
Aa_maderoi_ERR7622624_R1_fastp.fq
Aa_maderoi_ERR7622624_R2_fastp.fq
Anacamptis_pyramidalis_ERR7622662_R1_fastp.fq
Anacamptis_pyramidalis_ERR7622662_R2_fastp.fq
Apostasia_nuda_ERR5006156_R1_fastp.fq
Apostasia_nuda_ERR5006156_R2_fastp.fq
Coryanthes_hunteriana_SRR26876810_R1_fastp.fq
Coryanthes_hunteriana_SRR26876810_R2_fastp.fq
Cycnoches_amparoanum_SRR26876932_R1_fastp.fq
Cycnoches_amparoanum_SRR26876932_R2_fastp.fq
Dioscorea_caucasica_ERR4180009_R1_fastp.fq
Dioscorea_caucasica_ERR4180009_R2_fastp.fq
Paphiopedilum_callosum_SRR6127593_R1_fastp.fq
Paphiopedilum_callosum_SRR6127593_R2_fastp.fq
Vanilla_hartii_SRR26876866_R1_fastp.fq
Vanilla_hartii_SRR26876866_R2_fastp.fq
```

El recorte se ejecutó con el siguiente bucle:

```bash
while read F; do fastp -i ${F}_1.fastq.gz -I ${F}_2.fastq.gz -o ../trimmed/${F}_R1_fastp.fq -O ../trimmed/${F}_R2_fastp.fq -q 30 -w 8 -R ${F}.report -h ../trimmed/${F}.webpage; done < namelist.temp
```

Este comando:

- lee los nombres de muestra desde `namelist.temp`
- recorta lecturas pareadas
- escribe las lecturas recortadas en `../trimmed/`
- mantiene lecturas con umbral de calidad `q=30`
- genera un informe HTML por muestra

---

## 5. Ensamblaje de captura dirigida y recuperación de secuencias

### 5.1 Ensamblaje de loci con `HybPiper`

Activa primero el entorno conda de HybPiper:

```bash
conda activate hybpiper
```

Después ejecuta el bucle de ensamblaje:

```bash
while read F; do hybpiper assemble --readfiles ${F}_R1_fastp.fq ${F}_R2_fastp.fq --targetfile_dna /media/sf_VM_phylogenomics_Workshop_2026/data/HybPiper/RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta --cpu 2 --cov_cutoff 8 --bwa --timeout_assemble 80 --timeout_extract_contigs 100 --hybpiper_output ~/Documents/Phylo_work_2026/HybPiper --force_overwrite --prefix $F; done < namelist.temp
```

### 5.2 Resumen de recuperación de loci

Tras el ensamblaje, resume las estadísticas de recuperación con:

```bash
hybpiper stats --targetfile_dna RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta gene namelist.temp
```

Este paso ayuda a evaluar el éxito del enriquecimiento y del ensamblaje entre muestras y loci.

### 5.3 Recuperación de secuencias ensambladas

Recupera las secuencias de ADN para todas las muestras y loci con:

```bash
hybpiper retrieve_sequences dna --targetfile_dna RefGenes/OrchidsAspargales_Angio353_ITS_matK_rbcL_Orchs_targetfile.fasta --sample_names namelist.temp
```

Esto produce archivos FASTA por gen para los alineamientos posteriores.

---

## 6. Alineamiento por locus y filtrado

### 6.1 Alineamiento de loci con `MAFFT`

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

---

## 7. Inferencia de árboles génicos y postprocesamiento

### 7.1 Inferencia de árboles génicos con `RAxML-NG`

Crea primero una carpeta de salida:

```bash
mkdir -p ../raxml_gene_trees
```

Después ejecuta todos los alineamientos filtrados en paralelo:

```bash
ls *.fasta | parallel -j 2 'raxml-ng --all --msa {} --model GTR+G --bs-trees 100 --threads 2 --seed 4321 --prefix ../raxml_gene_trees/{.}'
```

Este comando:

- ejecuta un análisis de RAxML-NG por locus
- infiere un árbol de máxima verosimilitud bajo `GTR+G`
- realiza 100 réplicas de bootstrap
- produce árboles por locus con soportes de bootstrap

### 7.2 Organización de las salidas de RAxML-NG

Crea carpetas para ordenar las salidas:

```bash
mkdir -p output_byproduct MLtrees_support
```

Después mueve los archivos:

```bash
mv *raxml.* output_byproduct/
mv *.support MLtrees_support/
```

### 7.3 Renombrar y concatenar árboles con soporte

Dentro de la carpeta que contiene los árboles con soporte, renómbralos de `*.support` a `*.support.tre`:

```bash
rename 's/\.support$/.support.tre/' *.support
```

Después concaténalos en un único archivo:

```bash
cat *.tre > Orchidaceae_allSub_63g_MLtree.tre
```

### 7.4 Colapsar ramas con soporte débil

Usa Newick Utilities para colapsar ramas con soporte bootstrap menor o igual a 20:

```bash
nw_ed Orchidaceae_allSub_63g_MLtree.tre 'i & b<=20' o > Orchidaceae_allSub_63g_MLtree_20LBScoll.tre
```

Esto produce un conjunto filtrado de árboles génicos adecuado para análisis posteriores.

## 8. Inferencia del árbol de especies con ASTRAL

Una vez generados los árboles génicos y colapsadas las ramas débilmente soportadas, es posible inferir un **árbol de especies** usando **ASTRAL**, un método basado en cuartetos que estima el árbol de especies a partir de un conjunto de árboles génicos.

### 8.1 Preparar el conjunto de árboles génicos

En este taller, el archivo de entrada para ASTRAL es el conjunto de árboles génicos concatenados y filtrados:

```bash

Orchidaceae_allSub_63g_MLtree_20LBScoll.tre

```

Este archivo contiene todos los árboles génicos en formato Newick, con ramas de soporte bajo ya colapsadas.

### 8.2 Ejecutar ASTRAL

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

### 8.3 Interpretación general del resultado

El árbol resultante:

- representa una **hipótesis de relaciones entre especies**

- resume la señal compartida entre múltiples árboles génicos

- puede ayudar a identificar conflicto genealógico y discordancia entre loci

### 8.4 Salida esperada

El archivo principal será:

```bash

ASTRAL/Orchidaceae_allSub_63g_ASTRAL.tre

```

Este árbol puede visualizarse posteriormente en programas como **FigTree**, **iTOL** o **Dendroscope**.

### 8.5 Versión opcional con archivo de log

Si quieres guardar también el log de ejecución, puedes correr:

```bash

java -jar astral.jar \\

  -i Orchidaceae_allSub_63g_MLtree_20LBScoll.tre \\

  -o ASTRAL/Orchidaceae_allSub_63g_ASTRAL.tre \\

  > ASTRAL/astral.log

```

### 8.6 ¿Por qué usar ASTRAL?

ASTRAL es especialmente útil cuando:

- distintos genes apoyan historias evolutivas diferentes

- existe discordancia entre árboles génicos

- se desea inferir un árbol de especies a partir de múltiples loci independientes

En este contexto, ASTRAL complementa la inferencia de árboles génicos por **RAxML-NG** y permite pasar de una colección de gene trees a una hipótesis de árbol de especies.

---

## 9. Notas para estudiantes

- Comprueba siempre que estás en el directorio correcto con `pwd`.
- Usa `ls` con frecuencia para seguir la pista de tus archivos.
- Crea las carpetas de salida antes de redirigir archivos hacia ellas.
- Guarda los entornos y programas en el sistema de archivos local de Ubuntu y no dentro de carpetas compartidas de VirtualBox.
- Para trabajos largos, guarda los comandos en scripts de shell para poder relanzarlos fácilmente.
- Para resolver problemas, inspecciona la salida de `less`, `cat` y `grep` antes de borrar archivos.

---

## 10. Agradecimientos

La introducción breve a Bash usada en este taller fue adaptada del tutorial de Bash empleado en el repositorio del taller de ONT del mismo autor y reformulada aquí para una práctica de filogenómica.
